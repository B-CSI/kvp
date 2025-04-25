'''
KVP: Multiscale kurtosis phase picking
======================================

Kurtosis-Value-Picker (**KVP**) is a **seismic phase picker based on kurtosis** 
that produces picks with spectral information. It **filters input traces using 
a family of Ricker wavelet frames** to apply the wavelet transform (FIR filter),  
**achieving narrowband picking resolution**. The full algorithm is described 
in detail by the corresponding publication.

.. admonition:: Citation
    :class: important
    
    https://doi.org/10.1093/gji/ggaf136

.. currentmodule:: kvp.api

B-CSI implementation and Python package
---------------------------------------

Our implementation **revolves around the** :py:class:`KVP` **class**, which 
aims to provide a **simple and intuitive interface** to the algorithm. 
Instances of this class **store all picking parameters to run the full 
algorithm** on any data fed to them.

Picking **results are stored on instances of the** :py:class:`KVPOutput` 
**class**. This class provides access to **picked phases, available as POSIX 
timestamps or time in seconds from the start of the input data**. Optionally, 
this class **can also expose copies of both filtered and characteristic 
function (CF) traces**. This should be done carefully, as running the algorithm 
iterativelly over many data may quickly fill all available memory on your 
workstation.

.. admonition:: Memory optimizations
    :class: note
    
    The package tries to be smart about memory usage. :py:class:`KVP` instances 
    will try to reuse their already allocated memory during a previous run and 
    will only reallocate if necessary. This is all handled internally. In 
    particular, the performance boost when all inputs have the exact same 
    length (e.g. DAS) is massive.

'''


__version__ = '0.1.0'

DOI = '10.1093/gji/ggaf136'

BIBTEX = ''


from kvp.api import KVP, centralfreqs
