F2PY INSCTRUCTIONS:
==== =============
1) Create dynamic library from a fortran file:

    $ f2py -c cos_corr.f90 -m fast_corr
 
    (If the intel compiler is available, then 
    $ f2py -c cos_corr.f90 --fcompiler=intelem -m fast_corr  
    )

This will create a new file 'fort_lib.so' wich is callable in any python script

2) run python script normally:

   $ python get_lp_melt.py

   Output files are 'mean_Rnext_corr.mide' and 'mean_lp.mide'

   mean_Rnext_corr.mide gives the spatial orientational correlation along the chain. The first column is the distance 
   in number of bonds. The second is <cos(theta_s)>, where theta_s is the angle between bond i and bond i+s. The 
   < > average is performed over all bonds i

   The file mean_lp.mide gives the persistence length of the monomer in Lennard-Jonnes units
   l_p = - a / log < cos ( theta_1 ) >  
