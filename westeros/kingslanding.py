import numpy as np
from scipy.integrate import quad

def flux_converter(e_ref, e_min, e_max, index, flux=None, eflux=None, dnde=None, e2dnde=None):
    
    if flux!=None and eflux==None and dnde==None and e2dnde==None:
        N0 = flux/quad(lambda x: (x/e_ref)**(-index), e_min, e_max)[0]
        eflux = N0*quad(lambda x: x*(x/e_ref)**(-index), e_min, e_max)[0]
        dnde = N0
        e2dnde = e_ref**2*N0
        return flux, eflux, dnde, e2dnde
    if flux==None and eflux!=None and dnde==None and e2dnde==None:
        N0 = eflux/quad(lambda x: x*(x/e_ref)**(-index), e_min, e_max)[0]
        flux = N0*quad(lambda x: (x/e_ref)**(-index), e_min, e_max)[0]
        dnde = N0
        e2dnde = e_ref**2*N0
        return flux, eflux, dnde, e2dnde
    if flux==None and eflux==None and dnde!=None and e2dnde==None:
        N0 = dnde
        e2dnde = e_ref**2*N0
        flux = N0*quad(lambda x: (x/e_ref)**(-index), e_min, e_max)[0]
        eflux = N0*quad(lambda x: x*(x/e_ref)**(-index), e_min, e_max)[0]
        return flux, eflux, dnde, e2dnde
    if flux==None and eflux==None and dnde==None and e2dnde!=None:
        N0 = e2dnde/e_ref**2
        dnde = e2dnde/e_ref**2
        flux = N0*quad(lambda x: (x/e_ref)**(-index), e_min, e_max)[0]
        eflux = N0*quad(lambda x: x*(x/e_ref)**(-index), e_min, e_max)[0]
        return flux, eflux, dnde, e2dnde
    return 'Please give one and only one value in flux, eflux, dnde, and e2dnde.'
