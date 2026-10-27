'''
KVP: Multiscale kurtosis phase picking
======================================

Kurtosis-Value-Picker (**KVP**) is a **seismic phase picker based on kurtosis** 
that produces picks with spectral information. It **filters input traces using 
a family of Ricker or Morlet wavelet frames**, applying the wavelet transform 
(FIR filter) and **achieving narrowband picking resolution**. The full 
algorithm is detailed in the corresponding publication.

.. admonition:: Citation
    :class: important
    
    Hugo Latorre, Sergi Ventosa, Arantza Ugalde, Antonio Villaseñor, 
    Rafael Bartolomé, César R Ranero, KVP: a multiscale kurtosis approach for 
    seismic phase picking, Geophysical Journal International, Volume 241, 
    Issue 3, June 2025, Pages 1923–1935, https://doi.org/10.1093/gji/ggaf136

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
iterativelly over many data may quickly fill all available memory on the system.

.. admonition:: Memory optimizations
    :class: note
    
    The package tries to be smart about memory usage. :py:class:`KVP` instances 
    will try to reuse their already allocated memory during a previous run and 
    will only reallocate if necessary. This is all handled internally. In 
    particular, the performance boost when all inputs have the exact same 
    length (e.g. DAS) is massive.

'''


__version__ = '0.2.0'

DOI = '10.1093/gji/ggaf136'

BIBTEX = '@article{KVP,\n\
    author = {Latorre, Hugo and Ventosa, Sergi and Ugalde, Arantza and Villaseñor, Antonio and Bartolomé, Rafael and Ranero, César R},\n\
    title = {KVP: a multiscale kurtosis approach for seismic phase picking},\n\
    journal = {Geophysical Journal International},\n\
    volume = {241},\n\
    number = {3},\n\
    pages = {1923-1935},\n\
    year = {2025},\n\
    month = {04},\n\
    issn = {1365-246X},\n\
    doi = {10.1093/gji/ggaf136},\n\
    url = {https://doi.org/10.1093/gji/ggaf136},\n\
    eprint = {https://academic.oup.com/gji/article-pdf/241/3/1923/62907987/ggaf136.pdf},\n\
}'


from kvp.api import KVP, centralfreqs, W0
