import numpy as np
from scipy.integrate import quad

def angsep(lon1, lat1, lon2, lat2):
    
    '''
    angsep(lon1, lat1, lon2, lat2)
        
    Angular separation between two points (angles in degrees)
    
      The angular separation is calculated using the Vincenty formula [1]_,
      which is slighly more complex and computationally expensive than
      some alternatives, but is stable at at all distances, including the
      poles and antipodes.
      [1] http://en.wikipedia.org/wiki/Great-circle_distance
    '''
    
    deg2rad = np.pi/180
    lon1 = lon1*deg2rad
    lat1 = lat1*deg2rad
    lon2 = lon2*deg2rad
    lat2 = lat2*deg2rad

    sdlon = np.sin(lon2 - lon1)
    cdlon = np.cos(lon2 - lon1)
    slat1 = np.sin(lat1)
    slat2 = np.sin(lat2)
    clat1 = np.cos(lat1)
    clat2 = np.cos(lat2)

    num1 = clat2 * sdlon
    num2 = clat1 * slat2 - slat1 * clat2 * cdlon
    denom = slat1 * slat2 + clat1 * clat2 * cdlon

    return np.arctan2(np.sqrt(num1*num1 + num2*num2), denom)/deg2rad

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
