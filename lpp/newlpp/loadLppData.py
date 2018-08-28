#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 24 11:58:23 2018

Classes to read in data files containing either the 
light curve containing the transit (Timeseries) or
the information about the mapping (MapInfo)
These are used by lpp_transform

@author: smullally
"""

import scipy.io as spio
from astropy.io import fits
import requests

class TCE(object):
    
    def __init__(self, starid, ext=1,ddir=""):
        """
        starid is integer id, usually kicid
        """
        self.starid=starid
        self.filename = "%skplr%09u-20160128150956_dvt.fits" % (ddir,int(starid))
        self.ext=ext
        
        print self.filename        
        
    def readDV(self):

        try:
            hdu=fits.open(self.filename)
        except IOError:
            print "Filename not found."
            
        ext=self.ext
        
        self.time=hdu[ext].data['TIME']
        self.phase=hdu[ext].data['PHASE']
        self.flux=hdu[ext].data['LC_DETREND']
        self.period=hdu[ext].header['TPERIOD']
        self.tzero=hdu[ext].header['TEPOCH']
        self.dur=hdu[ext].header['TDUR']
        self.depth=hdu[ext].header['TDEPTH']
        self.mes=hdu[ext].header['MAXMES']
        
        hdu.close()
    
    def tceAPI(self):
        """
        Get all the data via the MAST API
        """
        url = "https://mastdev.stsci.edu"
        loc = 'api/v0.1/dvdata/%u/table/?tce=%u' % (self.starid,self.ext)
        getRequest= url + loc
        r=requests.get(url=getRequest)
        tce=r.json()
        
        self.time=getColumn(tce,'TIME')
        self.phase=getColumn(tce,'PHASE')
        self.flux=getColumn(tce,'LC_DETREND')
        self.period=
    
    def getColumn(self,tce,colname):
        data=np.array(map( lambda x : tce['data'][x][colname],\
                          np.arange(0,len(tce2['data']),1)))
        return data

class MapInfo(object):
    
    def __init__(self,filename):
        
        self.filename=filename
        
        self.readMatlabBlob(filename)

    
    def readMatlabBlob(self,filename):
        """
        read in matlab blob
        Using the DV trained one.
        """      

        mat=spio.loadmat(filename,matlab_compatible=True)
        
        #Pull out the information we need.
        
        self.n_dim = mat['mapInfoDV']['nDim'][0][0][0][0]
        self.Ymap = mat['mapInfoDV']['Ymap'][0][0][0][0]
        self.YmapMapping = self.Ymap['mapping']
        self.YmapMean = self.YmapMapping['mean'][0][0][0]
        self.YmapM = self.YmapMapping['M'][0][0]
        self.YmapMapped = self.Ymap['mapped']
        self.knn=mat['mapInfoDV']['knn'][0][0][0][0]
        self.knnGood=mat['mapInfoDV']['knnGood'][0][0][:,0]
        self.mappedPeriods=mat['mapInfoDV']['periods'][0][0][0]
        self.mappedMes=mat['mapInfoDV']['mes'][0][0][0]
        self.nPsample=mat['mapInfoDV']['nPsample'][0][0][0][0]  #number to sample
        self.nPercentil=mat['mapInfoDV']['npercentilTM'][0][0][0][0]
        self.dymeans=mat['mapInfoDV']['dymean'][0][0][0]
        self.ntrfr= 2.0
        self.npts=80.0
