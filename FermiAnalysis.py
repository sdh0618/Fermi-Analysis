#!/usr/bin/env python
# coding: utf-8

import numpy as np
import gt_apps

from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits

import make4FGLxml

class FermiAnalysis:



    def __init__(self, roi, path, evfile, scfile,
                 zmax=90, tmin=239557417, tmax=618049985, 
                 emin=500, emax=500000, ra=None, dec=None, lon=None, lat=None,
                 width=40, binsz=0.2, enumbins=None, ebinsperdec=8):
        
        self.roi = '{}_{}deg_{}MeV_{}GeV_{}_{}'.format(roi, width, emin, emax//1000, tmin, tmax)
        self.path = path
        self.evfile = evfile
        self.scfile = scfile
        
        self.zmax = zmax
        self.tmin = tmin
        self.tmax = tmax
        self.emin = emin
        self.emax = emax
        
        if (ra != None) & (dec != None) & (lat == None) & (lon == None):
            self.ra = ra
            self.dec = dec
            self.lon, self.lat = [float(x) for x in SkyCoord(ra*u.degree, dec*u.degree, frame='icrs').galactic.to_string().split()]
        elif (ra == None) & (dec == None) & (lat != None) & (lon != None):
            self.lon = lon
            self.lat = lat
            self.ra, self.dec = [float(x) for x in SkyCoord(lon*u.degree, lat*u.degree, frame='galactic').icrs.to_string().split()]
        else:
            print('Please provide the coordinates of the ROI center in ra/dec or lon/lat')
            pass
        
        self.width=width
        self.binsz=binsz
        self.enumbins=enumbins or round(ebinsperdec*np.log10(emax/emin))
            
        
        self.evfile_filterd = self.path + self.roi + '_filterd.fits'
        self.evfile_filterd_gti = self.path + self.roi + '_filterd_gti.fits'
        self.cmap = self.path + self.roi + '_cmap.fits'
        self.ltcube = self.path + self.roi + '_ltcube.fits'
        self.bexpmap = self.path + self.roi + '_bexpmap.fits'
        self.srcmdl = self.path + self.roi + '_srcmdl.xml'
        self.srcmaps = self.path + self.roi + '_srcmaps.fits'
    
    def gtselect(self, evclass=128, evtype=3, 
                 ra=None, dec=None, rad=None, emin=None, emax=None, zmax=None, tmin=None, tmax=None, 
                 infile=None, outfile=None,
                 chatter=4, clobber='yes', debug='no'):
                
        gt_apps.filter['evclass'] = evclass #Source=128, Clean=256, ultraclean=512, UltraCleanVeto=1024
        gt_apps.filter['evtype'] = evtype # (front + back)=3, front=1, back=2 PSF0=4, PSF1=8, PSF2=16, PSF3=32 (PSF1+PSF2+PSF3)=8+16+32=56
        gt_apps.filter['ra'] = ra or self.ra
        gt_apps.filter['dec'] = dec or self.dec
        gt_apps.filter['rad'] = rad or np.sqrt(2)*self.width/2
        gt_apps.filter['emin'] = emin or self.emin
        gt_apps.filter['emax'] = emax or self.emax
        gt_apps.filter['zmax'] = zmax or self.zmax
        gt_apps.filter['tmin'] = tmin or self.tmin
        gt_apps.filter['tmax'] = tmax or self.tmax
        gt_apps.filter['infile'] = infile or self.evfile
        gt_apps.filter['outfile'] = outfile or self.evfile_filterd
        gt_apps.filter['chatter'] = chatter
        gt_apps.filter['clobber'] = clobber
        gt_apps.filter['debug'] = debug
        gt_apps.filter.run()
        
    def gtmktime(self, filter='(DATA_QUAL>0)&&(LAT_CONFIG==1)', roicut='no',
                 scfile=None, evfile=None, outfile=None,
                 overwrite='no', chatter=4, clobber='yes', debug='no'):
                
        gt_apps.maketime['scfile'] = scfile or self.scfile
        gt_apps.maketime['filter'] = '(DATA_QUAL>0)&&(LAT_CONFIG==1)'
        gt_apps.maketime['roicut'] = 'no'
        gt_apps.maketime['evfile'] = evfile or self.evfile_filterd
        gt_apps.maketime['outfile'] = outfile or self.evfile_filterd_gti
        gt_apps.maketime['overwrite'] = 'no'
        gt_apps.maketime['chatter'] = 4
        gt_apps.maketime['clobber'] = 'yes'
        gt_apps.maketime['debug'] = 'no'
        gt_apps.maketime.run()
    
    def gtbin(self, algorithm='CCUBE', evfile=None, outfile=None, scfile=None,
              coordsys='GAL', proj='CAR', axisrot=0,
              nxpix=None, nypix=None, binsz=None, xref=None, yref=None,
              ebinalg='LOG', emin=None, emax=None, enumbins=None,
              chatter=4, clobber='yes', debug='no'):
                
        gt_apps.evtbin['algorithm'] = 'CCUBE'
        gt_apps.evtbin['evfile'] = evfile or self.evfile_filterd_gti
        gt_apps.evtbin['outfile'] = outfile or self.cmap
        gt_apps.evtbin['scfile'] = scfile or self.scfile
        gt_apps.evtbin['nxpix'] = nxpix or int(self.width/self.binsz)
        gt_apps.evtbin['nypix'] = nypix or int(self.width/self.binsz)
        gt_apps.evtbin['binsz'] = binsz or self.binsz
        gt_apps.evtbin['coordsys'] = 'GAL'
        gt_apps.evtbin['xref'] = xref or self.lon
        gt_apps.evtbin['yref'] = yref or self.lat
        gt_apps.evtbin['axisrot'] = 0
        gt_apps.evtbin['proj'] = 'CAR'
        gt_apps.evtbin['ebinalg'] = 'LOG'
        gt_apps.evtbin['emin'] = emin or self.emin
        gt_apps.evtbin['emax'] = emax or self.emax
        gt_apps.evtbin['enumbins'] = enumbins or self.enumbins 
        gt_apps.evtbin['chatter'] = 4
        gt_apps.evtbin['clobber'] = 'yes'
        gt_apps.evtbin['debug'] = 'no'
        gt_apps.evtbin.run()
        
    def gtltcube(self, zmax=None, dcostheta=0.025, binsz=1,
                evfile=None, scfile=None, outfile=None,
                chatter=4, clobber='yes', debug='no'):
        
        self.ltcube = self.path + self.roi + '_ltcube.fits'
        
        gt_apps.expCube['evfile'] = evfile or self.evfile_filterd_gti
        gt_apps.expCube['scfile'] = scfile or self.scfile
        gt_apps.expCube['outfile'] = outfile or self.ltcube
        gt_apps.expCube['zmax'] = zmax or self.zmax
        gt_apps.expCube['dcostheta'] = 0.025
        gt_apps.expCube['binsz'] = 1
        gt_apps.expCube['chatter'] = 4
        gt_apps.expCube['clobber'] = 'yes'
        gt_apps.expCube['debug'] = 'no'
        gt_apps.expCube.run()
        
    def gtexpcube2(self, irfs='P8R3_SOURCE_V3', evtype=3,
                   infile=None, outfile=None, cmap='none',
                   nxpix=None, nypix=None, binsz=None, xref=None, yref=None,
                   coordsys='GAL', lon=None, lat=None, axisor=0, proj='CAR',
                   ebinalg='LOG', emin=None, emax=None, enumbins=None,
                   chatter=4, clobber='yes', debug='no'):
        
        gt_apps.gtexpcube2['irfs'] = 'P8R3_SOURCE_V3'
        gt_apps.gtexpcube2['evtype'] = 3
        gt_apps.gtexpcube2['infile'] = infile or self.ltcube
        gt_apps.gtexpcube2['outfile'] = outfile or self.bexpmap
        gt_apps.gtexpcube2['cmap'] = cmap
        gt_apps.gtexpcube2['nxpix'] = int(360/(binsz or self.binsz))
        gt_apps.gtexpcube2['nypix'] = int(180/(binsz or self.binsz))
        gt_apps.gtexpcube2['binsz'] = binsz or self.binsz
        gt_apps.gtexpcube2['coordsys'] = 'GAL'
        gt_apps.gtexpcube2['xref'] = lon or self.lon
        gt_apps.gtexpcube2['yref'] = lat or self.lat
        gt_apps.gtexpcube2['axisrot'] = 0
        gt_apps.gtexpcube2['proj'] = 'CAR'
        gt_apps.gtexpcube2['ebinalg'] = 'LOG'
        gt_apps.gtexpcube2['emin'] = emin or self.emin
        gt_apps.gtexpcube2['emax'] = emax or self.emax
        gt_apps.gtexpcube2['enumbins'] = enumbins or self.enumbins
        gt_apps.gtexpcube2['chatter'] = 4
        gt_apps.gtexpcube2['clobber'] = 'yes'
        gt_apps.gtexpcube2['debug'] = 'no'
        gt_apps.gtexpcube2.run()

    def makexml(self, sources, ft1=None, out=None, DRversion=3,
                GDfile="$(FERMI_DIR)/refdata/fermi/galdiffuse/gll_iem_v07.fits", GDname='gll_iem_v07',
                ISOfile="$(FERMI_DIR)/refdata/fermi/galdiffuse/iso_P8R3_SOURCE_V3_v1.txt",
                ISOname='iso_P8R3_SOURCE_V3_v1',
                normsOnly=False, extDir='', radLim=-1, maxRad=None, ExtraRad=10, sigFree=5,
                varFree=True, psForce=False, E2CAT=False, makeRegion=True, GIndexFree=False, wd='', oldNames=False):
        '''
        make the xml file
        arguments are:
        GDfile (str) -- optional, location and name of Galactic diffuse model to use
        GDname (str) -- optional, name of Galactic diffuse component to use in xml model
        ISOfile (str) -- optional, location and name of Isotropic diffuse template to use
        ISOname (str) -- optional, name of Isotropic diffuse component to use in xml model
        normsOnly (bool) -- optional, flag to only set normalizations parameters free
        extDir (str) -- optional, directory with extended source templates
        radLim (float) -- optional, radius in degrees from center of ROI beyond which source parameters are fixed
        maxRad (float) -- optional, absolute maximum radius beyond which sources are fixed, this may be necessary when doing binned analysis and a variable source beyond radLim would be set free but this source is beyond the boundaries of the square region used for the binned likelihood
        ExtraRad (float) -- optional, radius beyond ROI radius in event file out to which sources will be included with fixed parameters, defaul tof 10 is good for analyses starting around 100 MeV, but for higher energy fits this can be decreased
        sigFree (float) -- optional, average significance, using FITS catalog file, below which source parameters are fixed, even if within radLim.  This corresponds to TS_value if using the XML catalog file.
        varFree (float) -- optional, variability index above which source parameters are free, if beyond radLim and/or below sigFree only the normalization parameter is set free
        psForce (bool) -- optional, flag to force extended sources to be point sources
        E2CAT (bool) -- optional, flag to force use catalog names for extended sources (only matters if using catalog FITS file)
        makeRegion (bool) -- optional, flag to also generate ds9 region file
        GIndexFree (bool) -- optional, the Galactic diffuse is given a power-law spectral shape but the by default the index is frozen, setting this flag to True allows that to be free for additional freedom in diffuse fit
        oldNames (bool) -- optional, flag to use the naming convention from make1FGLxml.py and make2FGLxml.py with a leading underscore and no spaces
        '''        
        make4FGLxml.srcList(sources, ft1=ft1 or self.evfile_filterd_gti, out=out or self.srcmdl, DRversion=DRversion).makeModel(
            GDfile=GDfile, GDname=GDname, ISOfile=ISOfile, ISOname=ISOname, normsOnly=normsOnly, 
            extDir=extDir, radLim=radLim, maxRad=maxRad, ExtraRad=ExtraRad, 
            sigFree=sigFree, varFree=varFree, psForce=psForce, E2CAT=E2CAT, 
            makeRegion=makeRegion, GIndexFree=GIndexFree, wd=wd, oldNames=oldNames)
        
    def gtsrcmaps(self, evtype=3, irfs='P8R3_SOURCE_V3', ptsrc='yes',
                  scfile=None, expcube=None, cmap=None, srcmdl=None, bexpmap=None, outfile=None,
                  chatter=4, clobber='yes', debug='no'):
        gt_apps.srcMaps['scfile'] = scfile or self.scfile
        gt_apps.srcMaps['expcube'] = expcube or self.ltcube
        gt_apps.srcMaps['cmap'] = cmap or self.cmap
        gt_apps.srcMaps['srcmdl'] = srcmdl or self.srcmdl
        gt_apps.srcMaps['bexpmap'] = bexpmap or self.bexpmap
        gt_apps.srcMaps['outfile'] = outfile or self.srcmaps
        gt_apps.srcMaps['evtype'] = evtype
        gt_apps.srcMaps['irfs'] = 'P8R3_SOURCE_V3'
        gt_apps.srcMaps['ptsrc'] = ptsrc
        gt_apps.srcMaps['chatter'] = 4
        gt_apps.srcMaps['clobber'] = 'yes'
        gt_apps.srcMaps['debug'] = 'no'
        gt_apps.srcMaps.run()