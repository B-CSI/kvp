# -*- coding: utf-8 -*-
"""
Minimal interface to wrappers around S. Ventosa original functions.
"""



import gc
from time import time


import numpy as np


from kvp.callers import get_jump_traces, update_jump_traces



class KurtosisJumps:
    '''Kurtosis jumps with a sliding window algorithm.
    '''
    
    def __init__(self, cycles : float):
        
        self._PARAM_ = \
            {
            'cycles' : float(cycles)
            }
    
        self._ITEMS_ = {}
        
        jump = \
            {
            'n_processed' : 0,
            'n_reset' : 0,
            'time_total' : 0.0
            }
        
        trig = \
            {
            'n_processed' : 0,
            'n_reset' : -1,
            'time_total' : 0.0
            }
        
        self._STATS_ = \
            {
            'jump' : jump,
            'trig' : trig,
            }
    
    def calculate(self, cfs, clear=False):
        '''Calculate internal vectors with jump values.
        '''
        
        t_ini = time()
        
        if not self._ITEMS_ or not self.check_alloc(cfs) or clear:
            cf_cycles = cfs._PARAM_['cycles']
            self._new(cfs._ITEMS_, cf_cycles)
            self._STATS_['jump']['n_reset'] += 1
        else:
            self._update_trig(cfs._ITEMS_)
        
        t_fin = time()
        
        self._ITEMS_['jump_time'] = t_fin - t_ini
        
        jump_time = t_fin - t_ini
        self._ITEMS_['jump_time'] = jump_time
        
        self._STATS_['jump']['n_processed'] += 1
        self._STATS_['jump']['time_total'] += jump_time
    
    def check_alloc(self, cfs):
        '''Check if allocated space is enough. Reallocate if not.
        '''
        
        cf_items = cfs._ITEMS_['cf']
        cf_fs = cf_items['fs']
        cf_centralfreq = cf_items['centralfreq']
        cf_downsample = cf_items['downsample']
        
        if not np.all(cf_fs==self._ITEMS_['trig']['fs']):
            return False
        
        if not np.all(cf_centralfreq==self._ITEMS_['trig']['centralfreq']):
            return False
        
        if not np.all(cf_downsample==self._ITEMS_['trig']['downsample']):
            return False
        
        data_N = cf_items['data_N']
        
        if data_N > self._ITEMS_['trig']['ALLOC']['DATA_N']:
            return False
        
        return True
    
    def _new(self, kur_items, cf_cycles):
        '''Internal.
        '''
        
        if self._ITEMS_: # Force garbage collection
            self._ITEMS_.clear()
            gc.collect()
        
        cycles = self._PARAM_['cycles']
        cf_windows = kur_items['cf']['windows']
        
        trig_windows = np.array(cycles/cf_cycles*cf_windows, dtype=np.uint32)
        
        trig_N = kur_items['cf']['N']
        trig_data = tuple(np.zeros(n, dtype=np.double) for n in trig_N)
        
        trig_delay = np.array(np.round(cf_windows/2), dtype=np.int64)
        
        cf_fs = kur_items['cf']['fs']
        cf_downsample = kur_items['cf']['downsample'].astype(np.int64)
        cf_centralfreq = kur_items['cf']['centralfreq']
        
        data_N = kur_items['cf']['data_N']
        DATA_N = kur_items['cf']['ALLOC']['DATA_N']
        
        ALLOC = \
            {
            'DATA' : trig_data,
            'N' : trig_N,
            'DATA_N' : DATA_N,
            }
        
        trig_items = \
            {
            'data' : trig_data,
            'data_N' : data_N,
            'windows' : trig_windows,
            'N' : trig_N.copy(),
            'fs' : cf_fs,
            'downsample' : cf_downsample,
            'delay' : trig_delay,
            'centralfreq' : cf_centralfreq,
            'ALLOC' : ALLOC,
            }
        
        self._ITEMS_['trig'] = trig_items
        
        cf_struct = kur_items['kur']['cf_struct']
        
        jump_items = get_jump_traces(trig_items, cf_struct)
        self._ITEMS_['jump'] = jump_items
    
    def _update_trig(self, cf_items):
        '''Internal.
        '''
        
        data_N = cf_items['cf']['data_N']
        
        if data_N != self._ITEMS_['trig']['data_N']:
            
            trig_N = self._ITEMS_['trig']['N']
            cf_N = cf_items['cf']['N']
            DATA = self._ITEMS_['trig']['ALLOC']['DATA']
            
            trig_N[:] = cf_N[:]
            
            trig_data_new = tuple(data[:n] for data, n in zip(DATA, trig_N))
            
            self._ITEMS_['trig']['data'] = trig_data_new
         
        jump_struct = self._ITEMS_['jump']['jump_struct']
        cf_struct = cf_items['kur']['cf_struct']
        
        update_jump_traces(jump_struct, cf_struct)
    
    @staticmethod
    def _mask_triggers(triggered, mingap):
        '''Internal.
        '''
        
        mask = np.diff(triggered, prepend=-mingap-1) > mingap
        return triggered[mask]
    
    def find_triggers(self, jump : float, mingap : float):
        '''Find triggering samples.
        '''
        
        t_ini = time()
        
        pick_param = \
            {
            'jump' : jump,
            'mingap' : mingap,
            }
        
        trig_fs = self._ITEMS_['trig']['fs']
        trig_windows = self._ITEMS_['trig']['windows']
        trig_data = self._ITEMS_['trig']['data']
        
        pick_mingap = np.array(np.round(mingap*trig_fs), dtype=np.int64)
        np.maximum(trig_windows, pick_mingap,
                                           out=pick_mingap, casting='unsafe')
        
        triggers = (np.where(tr>jump)[0] for tr in trig_data)
        
        selected = (self._mask_triggers(trig_, gap)
                                for trig_, gap in zip(triggers, pick_mingap))
        
        
        
        trig_items = \
            {
            'param' : pick_param,
            'mingap' : pick_mingap,
            'windows' : trig_windows,
            'fs' : trig_fs,
            'selected' : tuple(selected),
            }
        
        t_fin = time()
        
        trig_time = t_fin - t_ini
        
        trig_items['trig_time'] = trig_time
        
        self._STATS_['trig']['n_processed'] += 1
        self._STATS_['trig']['time_total'] += trig_time
        
        return trig_items
   