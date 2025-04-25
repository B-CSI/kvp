# -*- coding: utf-8 -*-
"""
Public API.
"""



from datetime import datetime, timezone
from io import StringIO, SEEK_END


from kvp.filter import MexhatFilters
from kvp.cf import KurtosisFunctions
from kvp.trigger import KurtosisJumps
from kvp.picker import KurtosisPicker
from kvp.callers import allocate_capture
from kvp.lib import KVPError



UTC = timezone.utc
TSFORMAT = '%Y-%m-%dT%H:%M:%S.%f'



class KVP:
    '''Main class to run the KVP algorithm. Instances of this class are used to 
    **run the algorithm with the stored parameters**, as long as the input data 
    is suitable (i.e. sufficient sampling rate and number of samples). The most 
    basic **input requires a data array, its sampling rate and an optional 
    POSIX timestamp** to be used as reference, i.e. the start time of the data. 
    ObsPy :py:class:`Trace` instances can be simply be passed to the 
    :py:meth:`obspy` method, which will handle everything.
    
    Inputs are fed through the :py:meth:`run` method.
    
    .. admonition:: Parameter restrictions
        :class: warning
        
        The following restrictions apply when setting the parameters or trying 
        to run the algorithm:
            
            * Jump window length cannot be longer than CF window length, and 
              typically should be much smaller (e.g. 2 cycles vs 90 cycles).
            * Any input array must be long enough to fit all CF windows, and 
              thus the largest scale, i.e. lowest central frequency, determines 
              the minimum required duration. Users must account for this when 
              cutting their data for inputs, ideally leaving space on both 
              sides.
            * Sampling rates must be at least five times larger than the 
              highest central frequency, otherwise the convolving wavelet would 
              not really be a Ricker.
            
    
    :param freqmax: Central frequency of the highest scale
    :type freqmax: :py:class:`float`
    :param octaves: Number of octaves
    :type octaves: :py:class:`int`
    :param voices: Number of voices per octave
    :type voices: :py:class:`int`
    :param cf_cycles: Length of CF sliding windows, in cycles
    :type cf_cycles: :py:class:`float`
    :param jmp_cycles: Length of jump windows, in cycles
    :type jmp_cycles: :py:class:`float`
    :param jump: Kurtosis jump threshold for triggers
    :type jump: :py:class:`float`
    :param mingap: Minimum gap between triggers, in seconds
    :type mingap: :py:class:`float`
    :param nbands: Required number of triggered bands to preserve picks
    :type nbands: :py:class:`int`
    :param downsample: Flag to downsample lower scales, default is :code:`True`
    :type downsample: :py:class:`bool`, optional
    :param N_max: Maximum number of expected possible grouped picks, default is 
        :code:`1000000`
    :type N_max: :py:class:`int`, optional 
    '''
    
    def __init__(self,
                 freqmax : float,
                 octaves : int,
                 voices : int,
                 cf_cycles : float,
                 jmp_cycles : float,
                 jump: float,
                 mingap : float,
                 nbands : int,
                 downsample=True,
                 N_max=1000000,
                 ):
        
        if jmp_cycles > cf_cycles:
            errmsg = 'jump cycles cannot be larger than CF cycles'
            raise ValueError(errmsg)
        
        self._PARAM_ = \
            {
            'freqmax' : freqmax,
            'octaves' : octaves,
            'voices' : voices,
            'cf_cycles' : cf_cycles,
            'jmp_cycles' : jmp_cycles,
            'jump' : jump,
            'mingap' : mingap,
            'nbands' : nbands,
            'downsample' : downsample,
            }
        
        filt = MexhatFilters(freqmax, octaves, voices, downsample)
        cf = KurtosisFunctions(cf_cycles)
        trig = KurtosisJumps(jmp_cycles)
        pick = KurtosisPicker(nbands)
        
        self._ITEMS_ = \
            {
            'filt' : filt,
            'cf' : cf,
            'trig' : trig,
            'pick' : pick,
            }
        
        self.ALLOC_CAPTURE(N_max)
        
        s = octaves * voices - 1
        freqmin = freqmax / 2**(s/voices)
        min_duration = cf_cycles / freqmin
        
        self._UTIL_ = \
            {
            'min_duration' : min_duration,
            'n_runs' : 0,
            }
    
    def obspy(self, trace, **kwargs):
        '''Run the picking algorithm for an input ObsPy 
        :py:class:`~obspy.Trace`.
        
        :param trace: Input trace
        :type trace: :py:class:`obspy.Trace`
        :param kwargs: Keyword arguments, these are passed to :py:meth:`run`
        '''
        
        data = trace.data
        fs = trace.stats.sampling_rate
        posix_ref = trace.stats.starttime.timestamp
        
        netw = trace.stats.network
        sta = trace.stats.station
        chn = trace.stats.channel
        
        rec_id = f'{netw}.{sta}.{chn}'
        
        return self.run(data, fs, posix_ref, rec_id, **kwargs)
    
    def run(self, data, fs, posix_ref=0.0,
            rec_id='', wt_traces=False, cf_traces=False):
        '''Run the picking algorithm for an input data array.
        
        :param data: Input data array
        :type data: :py:class:`numpy.ndarray`
        :param fs: Sampling rate of input data
        :type fs: :py:class:`float`
        :param posix_ref: Reference POSIX timestamp (starting time of data), 
            defaults to 0.0 if not given (1970-01-01T00:00:00) 
        :type posix_ref: :py:class:`float`, optional
        :param rec_id: Text identifier for the input data, e.g.
            :code:`{network}.{station}.{channel_code}`
        :type rec_id: :py:class:`str`, optional
        :param wt_traces: Flag to include filtered traces on the output, 
            default is :code:`False`
        :type wt_traces: :py:class:`bool`, optional
        :param cf_traces: Flag to include CF traces on the output, default is 
            :code:`False`
        :type cf_traces: :py:class:`bool`, optional
        '''
        
        kvp_items = self._ITEMS_
        kvp_param = self._PARAM_
        
        filt = kvp_items['filt']
        cf = kvp_items['cf']
        trig = kvp_items['trig']
        pick = kvp_items['pick']
        
        input_duration = data.size / fs
        
        if input_duration < self._UTIL_['min_duration']:
            errmsg = 'input data segment is too short'
            raise KVPError(errmsg)
        
        filt.calculate(data, fs)
        cf.calculate(filt)
        trig.calculate(cf)
        
        jump = kvp_param['jump']
        mingap = kvp_param['mingap']
        
        trig_items = trig.find_triggers(jump, mingap)
        refined_items = pick.refine_picks(trig_items, cf, posix_ref)
        self.debug = refined_items
        selected_picks = pick.select_picks(refined_items)
        
        self._UTIL_['current_picks'] = selected_picks
        self._UTIL_['current_id'] = rec_id
        self._UTIL_['current_posix_ref'] = posix_ref
        self._UTIL_['n_runs'] += 1
        
        return KVPOutput(self, wt_traces, cf_traces)
    
    def ALLOC_CAPTURE(self, N_max):
        '''Explicitly allocate memory for the capture algorithm, which avoids 
        continuous reallocation during successive runs if the maximum number of 
        possible picks grows after each run. A large enough pre-allocation will 
        improve performance, if it lasts.
        
        The default value should be enough for any reasonable input data.
        
        :param N_max: Maximum number of expected possible grouped picks
        :type N_max: :py:class:`int`
        '''
        
        nbands_max = self._PARAM_['octaves'] * self._PARAM_['voices']
        
        grouped_items = allocate_capture(N_max, nbands_max)
        
        self._ITEMS_['pick']._ITEMS_['grouped'] = grouped_items
        self._ITEMS_['pick']._STATS_['capture']['n_reset'] += 1



class KVPOutput:
    '''Class to store KVP outputs. It contains :py:class:`KVPick` instances and 
    supports list-like indexing, i.e. :code:`output[idx]`, to access them. 
    It also supports iteration to yield picks, i.e. :py:func:`iter` function 
    and :code:`for pick in output: ...` iterator syntax.
    
    If filtered or CF traces are saved after a run, they can be retrieved from 
    the corresponding methods, :py:meth:`wt_traces` and :py:meth:`cf_traces`, 
    respectively. The returned tuples contain a series of two-tuples with the 
    data arrays and dictionaries containing sampling rate and starttime.
    '''
    
    def __init__(self, kvp_items, wt_traces, cf_traces):
        
        wt = ()
        cf = ()
        
        picks = kvp_items._UTIL_['current_picks']
        rec_id = kvp_items._UTIL_['current_id']
        posix_ref = kvp_items._UTIL_['current_posix_ref']
        
        if wt_traces:
            wt = kvp_items._ITEMS_['filt'].traces()
            wt = tuple((data.copy(),
                        {'sampling_rate':fs,'starttime':posix_ref}, fc)
                                                       for data, fs, fc in wt)
        
        if cf_traces:
            cf = kvp_items._ITEMS_['cf'].traces()
            cf = tuple((data.copy(),
                        {'sampling_rate':fs,'starttime':posix_ref+wl/2.0},wl)
                                                   for data, fs, wl in cf)
        
        self._ITEMS_ = \
            {
            'picks' : picks,
            'rec_id' : rec_id,
            'posix_ref' : posix_ref,
            'wt_traces' : wt,
            'cf_traces' : cf,
            }
    
    def __iter__(self):
        '''Iterate over picks.
        '''
        
        return iter(self.picks)
    
    def __len__(self):
        '''Len method.
        '''
        
        return len(self.picks)
    
    def __bool__(self):
        '''Bool method.
        '''
        
        return bool(self.picks)
    
    def __getitem__(self, name):
        '''List-like behaviour.
        '''
        
        return self.picks[name]
    
    def __str__(self):
        '''Print method.
        '''
        
        with StringIO() as out:
            self._write_to_buffer(out)
            strout = out.getvalue()
        
        return strout
    
    def _write_to_buffer(self, buff_obj, posix=True, isoformat=True,
                                                                     **kwargs):
        '''Write to a file-like buffer (e.g. file object, StringIO).
        '''
        
        if not self.picks:
            return
        
        rec_id = self._ITEMS_['rec_id']
        
        for pick in self.picks:
            
            ons = pick.onset(posix=posix, **kwargs)
            nb = pick.nb
            fl = pick.centralfreqs[0]
            fh = pick.centralfreqs[-1]
            
            if posix and isoformat:
                ons = datetime.fromtimestamp(ons, UTC)
                ons = ons.strftime(TSFORMAT)
            
            buff_obj.write(f'{rec_id}\t{ons}\t{nb}\t{fl:.2f}\t{fh:.2f}\n')
        
        pos = buff_obj.seek(0, SEEK_END) - 1
        buff_obj.truncate(pos) # Remove last newline
        buff_obj.seek(pos)
    
    def append_to_file(self, fout, **kwargs):
        '''Append picks to an open file, i.e. an object returned by the 
        :py:func:`open` function. Output columns are :code:`rec_id`, onset 
        time, number of bands, lowest band, and highest band, respectively.
        
        :param fout: Open file object
        :type fout: file object
        :param posix: Flag to get onsets as POSIX time, default is 
            :code:`True`. If not, get time in seconds from input data start
        :type posix: :py:class:`bool`, optional
        :param isoformat: Flag to write POSIX times as ISO datetime strings or 
            :py:class:`float` numbers, default is :code:`True`. Only relevant 
            if :code:`posix` is set to :code:`True`
        :type isoformat: :py:class:`bool`, optional
        '''
        
        self._write_to_buffer(fout, **kwargs)
    
    def wt_traces(self):
        '''List of filtered traces, if available.
        '''
        
        return self._ITEMS_['wt_traces']
    
    def cf_traces(self):
        '''List of CF traces, if available.
        '''
        
        return self._ITEMS_['cf_traces']
    
    @property
    def picks(self):
        '''List of picks.
        '''
        
        return self._ITEMS_['picks']



def centralfreqs(freqmax, octaves, voices):
    '''Utility function to get the central frequencies for a wavelet family, 
    given the relevant KVP parameters.
    
    :param freqmax: Central frequency of the highest scale
    :type freqmax: :py:class:`float`
    :param octaves: Number of octaves
    :type octaves: :py:class:`int`
    :param voices: Number of voices per octave
    :type voices: :py:class:`int`
    '''
    
    fc = tuple(round(freqmax/2**(s/voices),2) for s in range(octaves*voices))
    
    return fc
    
    