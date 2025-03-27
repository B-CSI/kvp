#
'''
Binary interfaces.
'''



from pathlib import Path
from platform import system


import fancytypes as ft



LIBNAME = \
    {
     'Linux' : 'libkvp-unix.so',
     'Darwin' : 'libkvp-maco.so',
     'Windows' : 'libkvp-bill.dll',
     }


class KVPError(Exception):
    '''Exception for KVP related errors.
    '''
    
    pass



@ft.cstruct(wrapped=False)
class t_WaveletFamily:
    wr : ft.pointer(ft.pointer(ft.real64))
    scale : ft.pointer(ft.real64)
    Ls : ft.pointer(ft.uint32)
    center : ft.pointer(ft.int32)
    downsample : ft.pointer(ft.uint32)
    Cpsi : ft.real64
    op1 : ft.real64
    Ns : ft.uint32



@ft.cstruct(wrapped=False)
class t_RWTvar:
    d : ft.pointer(ft.pointer(ft.real64))
    N : ft.pointer(ft.uint32)
    S : ft.uint32



@ft.cstruct(wrapped=False)
class t_Movarray:
    x : ft.pointer(ft.pointer(ft.real64))
    y : ft.pointer(ft.pointer(ft.real64))
    dims : ft.pointer(ft.uint32)
    L : ft.pointer(ft.uint32)
    N : ft.uint32



@ft.cstruct
class t_Singlepick:
    idx_ref : ft.int64
    idx_max : ft.int64
    idx_ons : ft.int64
    t_ons : ft.real64
    posix_ons : ft.real64
    dt_jump : ft.real64
    kv_jump : ft.real64
    grouped : ft.logical
    pad : ft.character[7] # Explicitly doing this prevents an annoying (an useless) NumPy alignment-related warning



# Tentative experimental layout for an improved, future version of KVP
@ft.cstruct
class t_Singlepick_exp:
    idx_ref : ft.uint32
    idx_ons : ft.uint32
    idx_end : ft.uint32
    idx_max : ft.uint32
    t_ons : ft.real64
    posix_ons : ft.real64
    kv_jump : ft.real64
    dt_jump : ft.real64
    kv_peak : ft.real64
    dt_peak : ft.real64



@ft.cstruct(wrapped=False)
class t_KVPicks:
    picks : ft.pointer(ft.pointer(t_Singlepick))
    idx_bands : ft.pointer(ft.int32)
    nbands : ft.pointer(ft.int32)
    nbands_max : ft.int32
    N_max : ft.int32



_lib_name = LIBNAME.get(system())

if not _lib_name:
    errmsg = f'KVP binary not available for {system()} system'
    raise KVPError(errmsg) from None

_lib_path = Path(__file__).resolve().parent / _lib_name

if not _lib_path.is_file():
    errmsg = f'missing KVP library {_lib_name} for {system()} system'
    raise KVPError(errmsg)


lib_kvp = ft.load(_lib_path)


lib_kvp.setscales0 = \
    ft.interface(
        ft.pointer(ft.real64),
        ft.real64,
        ft.uint32,
        ft.uint32,
        ft.real64,
        )

lib_kvp.setwaveletlength0 = \
    ft.interface(
        ft.pointer(ft.uint32),
        ft.pointer(ft.real64),
        ft.uint32,
        ft.uint32,
        returns=ft.uint32,
        )

lib_kvp.setsampling0 = \
    ft.interface(
        ft.pointer(ft.uint32),
        ft.uint32,
        ft.int32,
        ft.uint32,
        ft.uint32,
        ft.real64,
        ft.real64,
        ft.real64,
        )

lib_kvp.FillWaveletFamily = \
    ft.interface(
        ft.pointer(t_WaveletFamily),
        )

lib_kvp.real_1D_wavelet_dec = \
    ft.interface(
        ft.pointer(t_RWTvar),
        ft.pointer(ft.real64),
        ft.uint32,
        ft.pointer(t_WaveletFamily),
        returns=ft.int32,
        )

lib_kvp.movkurtosis = \
    ft.interface(
        ft.pointer(t_Movarray),
        )

lib_kvp.movjumps = \
    ft.interface(
        ft.pointer(t_Movarray),
        )

lib_kvp.captures = \
    ft.interface(
        ft.pointer(t_KVPicks),
        ft.pointer(ft.pointer(t_Singlepick)),
        ft.pointer(ft.int32),
        ft.pointer(ft.real64),
        ft.pointer(ft.pointer(ft.logical)),
        returns=ft.int32,
        )
