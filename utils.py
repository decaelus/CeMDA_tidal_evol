import numpy as np
from astropy.units import cds
import astropy.units as u

#These models come from Lopez & Fortney (2013 -- http://iopscience.iop.org/article/10.1088/0004-637X/792/1/1/meta)

def Lopez_core_radius(Mcore):
    MEarth = 1.*cds.Mgeo
    REarth = 1.*cds.Rgeo
    
    return ((Mcore/MEarth).value)**(0.25)*REarth

def Lopez_Renv(Mp, fenv=0.05, Fp=1., age=5e9*u.year):
    MEarth = 1.*cds.Mgeo
    fenv_baseline = 0.05
    Fp_baseline = 1.
    age_baseline = 5e9*u.year
    
    return 2.06*cds.Rgeo*((Mp/MEarth).value)**(-0.21)*\
            (fenv/fenv_baseline)**(0.59)*\
            (Fp/Fp_baseline)**(0.044)*\
            ((age/age_baseline).value)**(-0.18)

#Ignores the apparently small contribution from the radiative exterior
def Lopez_superearth_radius(Mp, fenv=0.05, Fp=1., age=5e9*u.year):
    Mcore = (1. - fenv)*Mp
#    Mcore = Mp
    Rcore = Lopez_core_radius(Mcore)
    
    Renv = Lopez_Renv(Mp, fenv=fenv, Fp=Fp, age=age)
    
    return (Rcore + Renv).value
