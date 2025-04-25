# -*- coding: utf-8 -*-
"""
Minimal interface to S. Ventosa wavelet library. Real transform with mexhat.
"""



import math
import gc
from time import time


import numpy as np


from kvp.callers import (create_wavelet_family, get_filtered_traces,
                         update_filtered_traces)


class MexhatFilters:
    '''FIR filters using Ricker wavelets (mexican hat).
    '''
    
    def __init__(self, freqmax : float, octaves : int, voices : int,
                 downsample=True):
        
        self._PARAM_ = \
            {
            'freqmax' : float(freqmax),
            'octaves' : int(octaves),
            'voices' : int(voices),
            'downsample' : bool(downsample),
            }
        
        self._ITEMS_ = {}
        
        self._STATS_ = \
            {
            'n_processed' : 0,
            'n_reset' : 0,
            'time_total' : 0.0
            }
    
    def calculate(self, data, fs, clear=False):
        '''Apply filter to data.
        '''
        
        t_ini = time()
        
        if fs/5 < self.freqmax:
            errmsg = 'invalid maximum frequency, must be a fifth of or ' \
                     'lower than sampling rate'
            raise ValueError(errmsg)
        
        data = data.astype(np.double)
        
        if not self._ITEMS_ or not self.check_alloc(data, fs) or clear:
            self._new(data, fs)
            self._STATS_['n_reset'] += 1
        else:
            self._update_traces(data)
        
        t_fin = time()
        
        filter_time = t_fin - t_ini
        self._ITEMS_['filter_time'] = filter_time
        
        self._STATS_['n_processed'] += 1
        self._STATS_['time_total'] += filter_time
    
    def check_alloc(self, data, fs):
        '''Check if allocated space is enough. Reallocate if not.
        '''
        
        if fs != self._ITEMS_['wf']['fs']:
            return False
        
        if data.size > self._ITEMS_['tr']['ALLOC']['DATA_N']:
            return False
        
        return True
    
    def _new(self, data, fs):
        '''Internal.
        '''
        
        if self._ITEMS_: # Force garbage collection
            self._ITEMS_.clear()
            gc.collect()
        
        freqmax = self._PARAM_['freqmax']
        octaves = self._PARAM_['octaves']
        voices = self._PARAM_['voices']
        continuous = 0 if self._PARAM_['downsample'] else 1
        
        wf_items = create_wavelet_family(octaves, voices, continuous, fs,
                                         freqmax)
        self._ITEMS_['wf'] = wf_items
        
        N = data.size
        
        tr_items = self._allocate_traces(N)
        self._ITEMS_['tr'] = tr_items
        
        wf_struct = wf_items['wf_struct']
        
        rwt_items = get_filtered_traces(tr_items, data, wf_struct)
        self._ITEMS_['rwt'] = rwt_items
    
    def _allocate_traces(self, N):
        '''Internal.
        '''
        
        fs = self._ITEMS_['wf']['fs']
        downsample = self._ITEMS_['wf']['downsample']
        
        tr_fs = np.array(fs/downsample, dtype=np.double)
        tr_N = np.array([math.ceil(N/k) for k in downsample], dtype=np.uint32)
        tr_data = tuple(np.zeros(n, dtype=np.double) for n in tr_N)
        
        ALLOC = \
            {
            'DATA' : tr_data,
            'N' : tr_N,
            'DATA_N' : N,
            }
        
        items = \
            {
            'fs' : tr_fs,
            'N' : tr_N.copy(),
            'data' : tuple(data[:n] for data, n in zip(tr_data, tr_N)),
            'data_N' : N,
            'ALLOC' : ALLOC,
            }
        
        return items
    
    def _update_traces(self, data):
        '''Internal.
        '''
        
        N = data.size
        
        if N != self._ITEMS_['tr']['data_N']:
            
            downsample = self._ITEMS_['wf']['downsample']
            tr_N = self._ITEMS_['tr']['N']
            DATA = self._ITEMS_['tr']['ALLOC']['DATA']
            
            tr_N_new = np.array([math.ceil(N/k) for k in downsample],
                                                            dtype=np.uint32)
            tr_N[:] = tr_N_new[:]
            tr_data_new = tuple(data[:n] for data, n in zip(DATA, tr_N))
            
            self._ITEMS_['tr']['data'] = tr_data_new
        
        rwt_struct = self._ITEMS_['rwt']['rwt_struct']
        wf_struct = self._ITEMS_['wf']['wf_struct']
        
        update_filtered_traces(rwt_struct, data, wf_struct)
        
    @property
    def freqmax(self):
        '''Central frequency of the highest scale.
        '''
        
        return self._PARAM_['freqmax']
    
    @property
    def octaves(self):
        '''Number of octaves.
        '''
        
        return self._PARAM_['octaves']
    
    @property
    def voices(self):
        '''Number of voices per octave.
        '''
        
        return self._PARAM_['voices']
    
    def traces(self):
        '''Get a tuple of (data, fs, fc) tuples.
        '''
        
        tr_data = self._ITEMS_['tr']['data']
        tr_fs = self._ITEMS_['tr']['fs']
        wf_fc = self._ITEMS_['wf']['centralfreq']
        
        output = ((data, fs, round(fc,2))
                  for data, fs, fc in zip(tr_data, tr_fs, wf_fc))
        
        return tuple(output)
    