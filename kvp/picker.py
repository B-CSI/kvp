# -*- coding: utf-8 -*-
"""
Refining and grouping of picks based on S. Ventosa MATLAB adventures.
"""



import gc
from time import time


import numpy as np


from kvp.callers import get_captures, update_captures
from kvp.lib import t_Singlepick



class KurtosisPicker:
    '''Refine and group a set of picks.
    '''
    
    def __init__(self, nbands : int):
        
        self._PARAM_ = \
            {
            'nbands' : nbands,
            }
        
        self._ITEMS_ = {}
        
        refine = \
            {
            'n_processed' : 0,
            'n_reset' : -1, # This has not allocated memory
            'time_total' : 0.0,
            }
        
        capture = \
            {
            'n_processed' : 0,
            'n_reset' : 0,
            'time_total' : 0.0,
            }
        
        select = \
            {
            'n_processed' : 0,
            'n_reset' : -1, # Same
            'time_total' : 0.0,
            }
        
        self._STATS_ = \
            {
            'refine' : refine,
            'capture' : capture,
            'select' : select,
            }
    
    @staticmethod
    def _refine_singlepicks(trigs, wlen, cf, cf_wlen, fs, posix_ref):
        '''Internal.
        '''
        
        picks = np.zeros(trigs.size, dtype=t_Singlepick)
        half = wlen // 2
        delay = int(round(cf_wlen/2))
        
        # This can be done in C with a bit of care with triggers near the end of the CFs
        # Doubt it makes much of a difference, the captures should be a way bigger bottleneck
        # They are, by a large margin
        
        for i in range(picks.size):
            
            idx_ref = trigs[i] - half
            
            # Discard this pick completely
            if idx_ref < 0:
                picks[i][0] = -1
                picks[i][7] = True
                continue
             
            segment = cf[idx_ref:idx_ref+wlen]
            
            idx_max = np.argmax(segment)
            idx_ons = np.argmin(segment)
            t_ons = (idx_ref+idx_ons+delay) / fs
            dt_jump = (idx_max-idx_ons) / fs
            kv_jump = segment[idx_max] - segment[idx_ons]
            
            picks[i]['idx_ref'] = idx_ref
            picks[i]['idx_max'] = idx_max
            picks[i]['idx_ons'] = idx_ons
            picks[i]['t_ons'] = t_ons
            picks[i]['posix_ons'] = posix_ref + t_ons
            picks[i]['dt_jump'] = dt_jump
            picks[i]['kv_jump'] = kv_jump
        
        return picks
    
    def _new(self, refined):
        '''Internal.
        '''
        
        if self._ITEMS_: # Force garbage collection
            self._ITEMS_.clear()
            gc.collect()
        
        grouped_items = get_captures(refined)
        
        self._ITEMS_['grouped'] = grouped_items
    
    def _update_capture(self, refined):
        '''Internal.
        '''
        
        grouped_items = self._ITEMS_['grouped']
        
        idx_p = update_captures(refined, grouped_items)
        
        if idx_p != grouped_items['N']: # This is expected to almost always happen
            
            ALLOC = grouped_items['ALLOC']
            
            grouped_items['pick_ptrs'] = ALLOC['PICK_PTRS'][:idx_p]
            grouped_items['idx_bands'] = ALLOC['IDX_BANDS'][:idx_p]
            grouped_items['nbands'] = ALLOC['NBANDS'][:idx_p]
            
            grouped_items['N'] = idx_p
    
    def check_alloc(self, refined):
        '''Check if allocated space is enough. Reallocate if not.
        '''
        
        nbands_max = refined['nbands_max']
        N_max = refined['N_max']
        
        if N_max > self._ITEMS_['grouped']['ALLOC']['NP']:
            return False
        
        # Less bands than allocated can reuse but is rare within the same script so I will not deal with it
        # Actually it is a terrible idea to try that, some pointer aritmetic relies on nbands_max
        # Tl;dr: do not change this
        if nbands_max != self._ITEMS_['grouped']['ALLOC']['NB']:
            return False
        
        return True
    
    def refine_picks(self, trig, cfs, posix_ref=0.0):
        '''Refine picks to find onset samples.
        '''
        
        t_ini = time()
        
        trig_windows = trig['windows']
        trig_fs = trig['fs']
        trig_selected = trig['selected']
        cfs_data = cfs._ITEMS_['cf']['data']
        cf_windows = cfs._ITEMS_['cf']['windows']
        trig_mingap = trig['mingap']
        
        iter_items = zip(trig_selected, trig_windows,
                                                 cfs_data, cf_windows, trig_fs)
        
        refined = (self._refine_singlepicks(trigs, wl, cf, cfwl, fs, posix_ref)
                                   for trigs, wl, cf, cfwl, fs in iter_items)
        refined = tuple(refined)
        
        npicks = np.array(tuple(picks.size for picks in refined),
                                                              dtype=np.int32)
        
        N_max = np.sum(npicks)
        
        capture_gaps = np.maximum(trig_windows/trig_fs, trig_mingap/trig_fs)
        capture_gaps /= 2.0
        
        centralfreq = cfs._ITEMS_['cf']['centralfreq']
        
        refined_items = \
            {
            'refined' : refined,
            'npicks' : npicks,
            'centralfreq' : centralfreq,
            'capture_gaps' : capture_gaps,
            'posix_ref' : posix_ref,
            'nbands_max' : npicks.size,
            'N_max' : N_max,
            }
        
        t_fin = time()
        
        refine_time = t_fin - t_ini
        
        refined_items['refine_time'] = refine_time
        
        self._STATS_['refine']['time_last'] = refine_time
        
        self._STATS_['refine']['time_total'] += refine_time
        self._STATS_['refine']['n_processed'] += 1
        
        return refined_items
    
    def select_picks(self, refined, clear=False):
        '''Groups and select picks based on number of bands triggered.
        '''
        
        t_ini = time()
        
        if not self._ITEMS_ or not self.check_alloc(refined) or clear:
            self._new(refined)
            self._STATS_['capture']['n_reset'] += 1
        else:
            self._update_capture(refined)
        
        t_mid = time()
        
        grouped_items = self._ITEMS_['grouped']
        
        pick_ptrs = grouped_items['pick_ptrs']
        idx_bands = grouped_items['idx_bands']
        nbands = grouped_items['nbands']
        centralfreq = np.flip(refined['centralfreq'])
        
        iter_items = zip(pick_ptrs, idx_bands, nbands)
        
        NB = self._PARAM_['nbands']
        
        selected_picks = tuple(KVPick(ptrs, centralfreq[idx_b[:nb]])
                                       for ptrs, idx_b, nb in iter_items
                                                                   if nb>=NB)
        
        t_fin = time()
        
        capture_time = t_mid - t_ini
        select_time = t_fin - t_mid
        
        self._STATS_['capture']['time_last'] = capture_time
        self._STATS_['select']['time_last'] = select_time
        
        self._STATS_['capture']['time_total'] += capture_time
        self._STATS_['capture']['n_processed'] += 1
        
        self._STATS_['select']['time_total'] += select_time
        self._STATS_['select']['n_processed'] += 1
        
        return selected_picks



class KVPick:
    '''Grouped KVP phase pick.
    '''
    
    def __init__(self, single_ptrs, bands):
        
        nb = bands.size
        picks_mem = t_Singlepick[nb]()
        
        for i in range(nb):
            picks_mem[i] = single_ptrs[i][0]
        
        self._nb_ = nb
        self._bands_ = LOCK_ARR(np.round(bands,2))
        self._PICKS_MEM_ = picks_mem
        self._picks_ = LOCK_ARR(np.ctypeslib.as_array(picks_mem,
                                                              t_Singlepick))
    
    def onset(self, method='med', posix=True):
        '''Get onset time.
        
        :param method: Method to calculate the onset, only median is 
            implemented for now (:code:`'med'`)
        :type method: :py:class:`str`
        :param posix: Flag to get onset as POSIX time, default is :code:`True`. 
            If not, get time in seconds from input data start
        :type posix: :py:class:`bool`, optional
        '''
        
        if method != 'med':
            errmsg = f'invalid onset picking method \'{method}\''
            raise ValueError(errmsg)
        
        idx = 4 if posix else 3
        
        times = tuple(sp[idx] for sp in self._picks_)
        
        onset = np.median(times)
        
        return onset
    
    @property
    def singlepicks(self):
        '''Individual singlepicks.
        '''
        
        return self._picks_
    
    @property
    def nb(self):
        '''Number of bands.
        '''
        
        return self._nb_
    
    @property
    def centralfreqs(self):
        '''Central frequencies.
        '''
        
        return self._bands_



def LOCK_ARR(arr):
    '''Lock array.
    '''
    
    arr.flags['WRITEABLE'] = False
    
    return arr