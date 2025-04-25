# -*- coding: utf-8 -*-
"""
Minimal interface to wrappers around S. Ventosa original functions.
"""



import gc
from time import time


import numpy as np


from kvp.callers import get_cf_traces, update_cf_traces



class KurtosisFunctions:
    '''Moving kurtosis with a stable sliding window algorithm.
    '''
    
    def __init__(self, cycles : float):
        
        self._PARAM_ = \
            {
            'cycles' : float(cycles)
            }
    
        self._ITEMS_ = {}
        
        self._STATS_ = \
            {
            'n_processed' : 0,
            'n_reset' : 0,
            'time_total' : 0.0
            }
    
    def calculate(self, filtered, clear=False):
        '''Calculate characteristic functions.
        '''
        
        t_ini = time()
        
        if not self._ITEMS_ or not self.check_alloc(filtered) or clear:
            self._new(filtered._ITEMS_)
            self._STATS_['n_reset'] += 1
        else:
            self._update_cf(filtered._ITEMS_)
        
        t_fin = time()
        
        cf_time = t_fin - t_ini
        self._ITEMS_['cf_time'] = cf_time
        
        self._STATS_['n_processed'] += 1
        self._STATS_['time_total'] += cf_time
    
    def check_alloc(self, filtered):
        '''Check if allocated space is enough. Reallocate if not.
        '''
        
        filt_items = filtered._ITEMS_
        wf_fs = filt_items['tr']['fs']
        wf_centralfreq = filt_items['wf']['centralfreq']
        wf_downsample = filt_items['wf']['downsample']
        
        if not np.all(wf_fs==self._ITEMS_['cf']['fs']):
            return False
        
        if not np.all(wf_centralfreq==self._ITEMS_['cf']['centralfreq']):
            return False
        
        if not np.all(wf_downsample==self._ITEMS_['cf']['downsample']):
            return False
        
        data_N = filtered._ITEMS_['tr']['data_N']
        
        if data_N > self._ITEMS_['cf']['ALLOC']['DATA_N']:
            return False
        
        return True
    
    def _new(self, filt_items):
        '''Internal.
        '''
        
        if self._ITEMS_: # Force garbage collection
            self._ITEMS_.clear()
            gc.collect()
        
        cycles = self._PARAM_['cycles']
        centralfreq = filt_items['wf']['centralfreq']
        tr_fs = filt_items['tr']['fs']
        tr_downsample = filt_items['wf']['downsample']
        
        cf_windows = np.array(cycles/centralfreq*tr_fs, dtype=np.uint32)
        
        cf_N = filt_items['tr']['N']
        cf_data = tuple(np.zeros(n, dtype=np.double) for n in cf_N)
        
        data_N = filt_items['tr']['data_N']
        DATA_N = filt_items['tr']['ALLOC']['DATA_N']
        
        ALLOC = \
            {
            'DATA' : cf_data,
            'N' : cf_N,
            'DATA_N' : DATA_N,
            }
        
        cf_items = \
            {
            'data' : tuple(data[:n] for data, n in zip(cf_data, cf_N)),
            'data_N' : data_N,
            'windows' : cf_windows,
            'N' : cf_N.copy(),
            'fs' : tr_fs,
            'centralfreq': centralfreq,
            'downsample' : tr_downsample,
            'ALLOC' : ALLOC,
            }
        
        self._ITEMS_['cf'] = cf_items
        
        rwt_struct = filt_items['rwt']['rwt_struct']
        
        kur_items = get_cf_traces(cf_items, rwt_struct)
        self._ITEMS_['kur'] = kur_items
    
    def _update_cf(self, filt_items):
        '''Internal
        '''
        
        data_N = filt_items['tr']['data_N']
        
        if data_N != self._ITEMS_['cf']['data_N']:
            
            cf_N = self._ITEMS_['cf']['N']
            tr_N = filt_items['tr']['N']
            DATA = self._ITEMS_['cf']['ALLOC']['DATA']
            
            cf_N[:] = tr_N[:]
            
            cf_data_new = tuple(data[:n] for data, n in zip(DATA, cf_N))
            
            self._ITEMS_['cf']['data'] = cf_data_new
         
        cf_struct = self._ITEMS_['kur']['cf_struct']
        rwt_struct = filt_items['rwt']['rwt_struct']
        
        update_cf_traces(cf_struct, rwt_struct)
    
    @property
    def cycles(self):
        '''Moving window lengths, in cycles.
        '''
        
        return self._PARAM_['cycles']
    
    def traces(self):
        '''Get tuple of (data, fs, wlen) tuples.
        '''
        
        cf_data = self._ITEMS_['cf']['data']
        cf_fs = self._ITEMS_['cf']['fs']
        cf_wlen = self._ITEMS_['cf']['windows']
        
        output = ((data, fs, wlen/fs)
                  for data, fs, wlen in zip(cf_data, cf_fs, cf_wlen))
        
        return tuple(output)
        