# -*- coding: utf-8 -*-
"""
Calls to binaries.
"""



import numpy as np
import fancytypes as ft


from kvp.lib import (t_WaveletFamily, t_RWTvar, t_Movarray,
                     t_Singlepick, t_KVPicks, lib_kvp)



# Mexhat fixed parameter
_OP1 = np.sqrt(2.0)
_A0 = 2.0
_B0 = 0.5
_N = 10000



def create_wavelet_family(octaves, voices, continuous, fs, freqmax):
    '''Internal.
    '''
    
    s0 = _OP1*fs / (2.0*np.pi*freqmax)
    Ns = octaves * voices
    
    scale = np.zeros(Ns, dtype=np.double)
    Ls = np.zeros(Ns, dtype=np.uint32)
    center = np.zeros(Ns, dtype=np.int32)
    downsample = np.zeros(Ns, dtype=np.uint32)
    
    scale_p = ft.real64.array(scale)
    Ls_p = ft.uint32.array(Ls)
    center_p = ft.int32.array(center)
    downsample_p = ft.uint32.array(downsample)
    
    lib_kvp.setscales0(scale_p, s0, Ns, voices, _A0)
    Nsmpl = lib_kvp.setwaveletlength0(Ls_p, scale_p, Ns, _N)
    lib_kvp.setsampling0(downsample_p, Ns, continuous, voices, octaves,
                        _B0, _A0, s0)
    
    frames = np.zeros(Nsmpl, dtype=np.double)
    
    frame_ptrs = (ft.pointer(ft.real64)*Ns)(ft.real64.array(frames))
    
    wf_struct = \
        t_WaveletFamily(
            frame_ptrs,
            scale_p,
            Ls_p,
            center_p,
            downsample_p,
            0.0,
            _OP1,
            Ns
            )
    
    lib_kvp.FillWaveletFamily(wf_struct)
    
    centralfreq = _OP1*fs / (2.0*np.pi*scale)
    
    items = \
        {
        'fs' : fs,
        'Ns' : Ns,
        'scale' : scale,
        'Ls' : Ls,
        'center' : center,
        'downsample' : downsample,
        'centralfreq' : centralfreq,
        'frames' : frames,
        'frame_ptrs' :frame_ptrs,
        'wf_struct' : wf_struct,
        }
    
    return items


def get_filtered_traces(tr_items, data_in, wf_struct):
    '''Internal.
    '''
    
    data_out = tr_items['ALLOC']['DATA']
    N_out = tr_items['N']
    
    data_ptrs = ft.ptrarray(data_out, ft.real64)
    
    rwt_struct = \
        t_RWTvar(
            data_ptrs,
            ft.uint32.array(N_out),
            N_out.size,
            )
    
    rwt_items = \
        {
        'data_ptrs' : data_ptrs,
        'rwt_struct' : rwt_struct,
        }
    
    N_in = data_in.size
    
    err = lib_kvp.real_1D_wavelet_dec(rwt_struct, ft.real64.array(data_in),
                                      N_in, wf_struct)
    
    return rwt_items


def update_filtered_traces(rwt_struct, data_in, wf_struct):
    '''Internal.
    '''
    
    N_in = data_in.size
    
    err = lib_kvp.real_1D_wavelet_dec(rwt_struct, ft.real64.array(data_in),
                                      N_in, wf_struct)


def get_cf_traces(cf_items, rwt_struct):
    '''Internal.
    '''
    
    data_out = cf_items['data']
    wlen = cf_items['windows']
    N_out = cf_items['N']
    
    data_ptrs = ft.ptrarray(data_out, ft.real64)
    
    kur_struct = \
        t_Movarray(
            rwt_struct.d,
            data_ptrs,
            ft.uint32.array(N_out),
            ft.uint32.array(wlen),
            N_out.size,
            )
    
    lib_kvp.movkurtosis(kur_struct)
    
    kur_items = \
        {
        'data_ptrs' : data_ptrs,
        'cf_struct' : kur_struct
        }
    
    return kur_items


def update_cf_traces(cf_struct, rwt_struct):
    '''Internal.
    '''
    
    cf_struct.x = rwt_struct.d
    
    lib_kvp.movkurtosis(cf_struct)


def get_jump_traces(trig_items, cf_struct):
    '''Internal.
    '''
    
    data_out = trig_items['data']
    wlen = trig_items['windows']
    N_out = trig_items['N']
    
    data_ptrs = ft.ptrarray(data_out, ft.real64)
    data_ptrs = ft.ptrarray(data_out, ft.real64)
    
    jump_struct = \
        t_Movarray(
            cf_struct.y,
            data_ptrs,
            ft.uint32.array(N_out),
            ft.uint32.array(wlen),
            N_out.size,
            )
    
    lib_kvp.movjumps(jump_struct)
    
    jump_items = \
        {
        'data_ptrs' : data_ptrs,
        'jump_struct' : jump_struct
        }
    
    return jump_items


def update_jump_traces(jump_struct, cf_struct):
    '''Internal.
    '''
    
    jump_struct.x = cf_struct.y
    
    lib_kvp.movjumps(jump_struct)


def get_captures(refined_items):
    '''I did it, and for no reason..
    '''
    
    refined = refined_items['refined']
    npicks = refined_items['npicks']
    capture_gaps = refined_items['capture_gaps']
    nbands_max = refined_items['nbands_max']
    N_max = refined_items['N_max']
    
    refined = tuple(reversed(refined))
    gaps = np.flip(capture_gaps).copy() # Needs to be a physical flip, not a view
    npicks_ = np.flip(npicks).copy() # Nasty out of bounds bug before I caught this (it used to be computed here from the reversed tuple)
    
    grouped_items = allocate_capture(N_max, nbands_max)
    
    capture_struct = grouped_items['capture_struct']
    
    refined_ptrs = ft.ptrarray(refined, t_Singlepick)
    
    mask = tuple(np.zeros(n, dtype=np.bool_) for n in npicks_) # np.bool_ is compatible with NumPy<2.0
    mask_ptrs = ft.ptrarray(mask, ft.logical)
    
    idx_p = lib_kvp.captures(capture_struct, refined_ptrs,
                             ft.int32.array(npicks_), ft.real64.array(gaps),
                             mask_ptrs)
    
    del refined_ptrs
    del mask_ptrs
    del mask
    
    if idx_p != grouped_items['N']: # This is expected to almost always happen
        
        ALLOC = grouped_items['ALLOC']
        
        grouped_items['pick_ptrs'] = ALLOC['PICK_PTRS'][:idx_p]
        grouped_items['idx_bands'] = ALLOC['IDX_BANDS'][:idx_p]
        grouped_items['nbands'] = ALLOC['NBANDS'][:idx_p]
        
        grouped_items['N'] = idx_p
    
    return grouped_items


def update_captures(refined_items, grouped_items):
    '''Internal.
    '''
    
    refined = refined_items['refined']
    npicks = refined_items['npicks']
    capture_gaps = refined_items['capture_gaps']
    
    capture_struct = grouped_items['capture_struct']
    
    refined = tuple(reversed(refined))
    gaps = np.flip(capture_gaps).copy()
    npicks_ = np.flip(npicks).copy()
    
    
    refined_ptrs = ft.ptrarray(refined, t_Singlepick)
    
    mask = tuple(np.zeros(n, dtype=np.bool_) for n in npicks_)
    mask_ptrs = ft.ptrarray(mask, ft.logical)
    
    idx_p = lib_kvp.captures(capture_struct, refined_ptrs,
                             ft.int32.array(npicks_), ft.real64.array(gaps),
                             mask_ptrs)
    
    del refined_ptrs
    del mask_ptrs
    del mask
    
    
    # Python fuction, not worth it even at low number of singlepicks
    # pick_ptrs = grouped_items['ALLOC']['PICK_PTRS']
    # idx_bands = grouped_items['ALLOC']['IDX_BANDS']
    # nbands = grouped_items['ALLOC']['NBANDS']
    
    # idx_p = capture_py(pick_ptrs, idx_bands, nbands, refined, npicks_, gaps)
    
    
    return idx_p # This needs to be updated outside, no grouped dict in here


def capture_py(pick_ptrs, idx_bands, nbands, refined, npicks, gaps):
    '''Python version of the capture algorithm. It might just be faster for 
    short traces (i.e. low number of picks). Actually, it is not. It stays here 
    because it describes the capture algorithm in detail.
    '''
    
    nbands_max = len(refined) # Number of frequency bands / scales
    
    idx_p = 0 # Pick index to store picks though the N_max axis
    idx_aux = np.zeros(nbands_max, np.int32) # Auxiliary index to keep track of searched picks
    
    pick_ptrs_aux = ft.pointer(t_Singlepick)[nbands_max]() # Auxiliary pointer array to group picks
    idx_bands_aux = np.zeros(nbands_max, dtype=np.int32) # Auxiliary band index array to group picks
    ptr_wiper = ft.pointer(t_Singlepick)[nbands_max]() # Auxiliary pointer array to wipe auxiliary into NULL references
    
    refined_ptrs = ft.ptrarray(refined, t_Singlepick)
    
    # Cascade from each band capturing picks
    for idx_b in range(nbands_max):
        
        ref_band = refined[idx_b] # "Point" to reference band
        gap = gaps[idx_b] # Use the corresponding gap time for captures
        
        # Search for reference picks within a band, skip if already grouped
        for idx_r in range(npicks[idx_b]):
            
            ref_pick = ref_band[idx_r] # Take a pick
            
            # Skip it if already grouped, relevant after completing the first band
            if ref_pick['grouped']:
                continue
            
            nb = 1 # Start a counter for number of captured picks (number of bands)
            
            # Band index and pick pointer into auxiliaries
            idx_bands_aux[0] = idx_b
            pick_ptrs_aux[0] = ft.cpointer(refined_ptrs[idx_b][idx_r])
            ref_pick['grouped'] = True # Flag as grouped
            
            ref_ons = ref_pick['t_ons'] # Take the reference onset
            
            # Cascade through remaining bands
            for idx_c in range(idx_b+1,nbands_max):
                
                cas_band = refined[idx_c] # "Point" to cascaded band
                search_from = idx_aux[idx_c] # Start from where the previous reference stopped (onsets are ordered)
                
                cap_idx = -1 # Just to catch bugs if this gets used somehow without setting, should never happen
                
                # Search through a band for captures
                for idx_s in range(search_from,npicks[idx_c]):
                    
                    s_pick = cas_band[idx_s] # Take a pick
                    
                    # Skip if grouped or before reference pick
                    if s_pick['grouped'] or ref_ons-s_pick['t_ons'] > gap:
                        continue
                    
                    # Stop search if after reference pick, setting the new start point for next references
                    if ref_ons-s_pick['t_ons'] < -gap:
                        idx_aux[idx_c] = idx_s
                        break
                    
                    # If none of those conditions is True then we got a potential capture
                    # Make sure to capture the closest pick to reference
                    
                    dt = abs(ref_ons-s_pick['t_ons']) # Starting time delta
                    cap_idx = idx_s # Store index of potential capture
                    
                    # Look for the minimum time delta
                    for idx_m in range(idx_s+1,npicks[idx_c]):
                        
                        m_pick = cas_band[idx_m] # Take the next pick
                        
                        # If grouped then it is not closer, it can be sometimes but requires undoing previous captures
                        # Not gonna bother now, it is a very niche issue
                        if m_pick['grouped']:
                            break
                        
                        # Stop capture if next pick is not closer
                        if abs(ref_ons-m_pick['t_ons']) > dt:
                            break
                        
                        # If the previous conditions are not met, then this new pick is closer
                        
                        dt = abs(ref_ons-m_pick['t_ons']) # Update time delta
                        cap_idx = idx_m # Update capture index
                    
                    # Store pick in auxiliaries and advance capture counter (number of band)
                    idx_bands_aux[nb] = idx_c
                    pick_ptrs_aux[nb] = \
                                    ft.cpointer(refined_ptrs[idx_c][cap_idx])
                    refined_ptrs[idx_c][cap_idx].grouped = True # Flag as grouped
                    nb += 1
                    
                    idx_aux[idx_c] = cap_idx # Set new start point for next references
                    
                    # Breaking at the loop level should be safe
                    # It can only happen if a pick is captured, otherwise it either stops iterating or breaks before
                    break
                
                # This is the cascading loop level, nothing happens here
            
            # This is the reference pick level, the grouped picks are dumped here
            nbands[idx_p] = nb
            idx_bands[idx_p][:] = idx_bands_aux[:]
            pick_ptrs[idx_p][:] = pick_ptrs_aux[:]
            
            # Wipe auxiliaries, number of bands gets restarted during the loop
            idx_bands_aux[:] = 0
            pick_ptrs_aux[:] = ptr_wiper[:]
            
            idx_p += 1 # Advance pick index
        
        # Reference band level, maybe set grouped to False so the grouping procedure can be called again
        # Or just make a boolean array on the spot for each band, that would have made more sense
        # Wipe the auxiliary index array so the next band can properly look for captures, sad times when I forgot this at first
        idx_aux[:] = 0
    
    # End of capture algorithm
    
    return idx_p


def allocate_capture(N_max, nbands_max):
    '''Manually allocate memory for captures. Number of singlepicks is not 
    really a predictable parameter, so having this allows an educated guess. 
    The performance gain could be noticeable.
    '''
    
    idx_bands = np.zeros((N_max,nbands_max), dtype=np.int32)
    nbands = np.zeros(N_max, dtype=np.int32)
    
    pick_ptrs = ft.pointer(t_Singlepick)[nbands_max][N_max if N_max else 1]()
    
    capture_struct = \
        t_KVPicks(
            pick_ptrs[0],
            ft.int32.array(idx_bands),
            ft.int32.array(nbands),
            nbands_max,
            N_max,
            )
    
    ALLOC = \
        {
        'PICK_PTRS' : pick_ptrs,
        'IDX_BANDS' : idx_bands,
        'NBANDS' : nbands,
        'NB' : nbands_max,
        'NP' : N_max
        }
    
    grouped_items = \
        {
        'pick_ptrs' : pick_ptrs[:N_max],
        'idx_bands' : idx_bands[:N_max],
        'nbands' : nbands[:N_max],
        'N' : N_max,
        'ALLOC' : ALLOC,
        'capture_struct' : capture_struct,
        }
    
    return grouped_items
