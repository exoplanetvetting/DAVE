#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  7 14:36:02 2018

@author: smullally

kwargs include ext=1 and header=True
"""

from astroquery.mast import Observations
from astropy.io import fits
from dave.fileio.AbstractMast import MastArchive

class TessAstroqueryArchive(MastArchive):
    
    def __init__(self, downloadDir="./mastDownload/"):
        
        self.path=downloadDir
        
        MastArchive.__init__(self, downloadDir, 'www.example.com')

        
    
    def downloadDVData(self):
        
        manifest=Observations.download_products(self.obsid, \
                              productSubGroupDescription=['DVT'], \
                              mrp_only=False, extension="fits",\
                              download_dir=self.path)

        self.dvFileLoc=manifest['Local Path'][0]

    
    def downloadLCData(self):
        manifest=Observations.download_products(self.obsid,\
                                            productSubGroupDescription=['LC'],\
                                            mrp_only=False,extension="fits",
                                            download_dir=self.path)
        
        self.lcFileLoc=manifest['Local Path'][0]

    def downloadTpfData(self):
        manifest=Observations.download_products(self.obsid,\
                                            productSubGroupDescription=['TP'],\
                                            mrp_only=False,extension="fits",
                                            download_dir=self.path)

        self.lcFileLoc=manifest['Local Path'][0]

    
    
#    def makeTessCut(self,size):
        
        

    def getTwoMinObservationId(self):
        """
        Using Astroquery find the two minute obsid.
        ticid = integer tic id
        sector = integer sector number
        """
        
        obs_query_string = "tess*-s%04i-%016i*" % (self.sector, self.ticid)
        obsTable = Observations.query_criteria(mission="TESS", obs_id=obs_query_string)
        
        self.obsid=obsTable['obsid']
        
        return obsTable
    
    def getDvt(self, ticid, sector, *args, **kwargs):
        """
        Combine the search and the download.
        """
        self.ticid=ticid
        self.sector=sector
        self.getTwoMinObservationId()
        self.downloadDVData()
        
        localUrl=self.dvFileLoc
        
        return self.parse(localUrl, *args, **kwargs)

    def getTpf(self, ticid, sector, *args, **kwargs):
        """
        Combine the search and the download.
        """
        self.ticid=ticid
        self.sector=sector
        self.getTwoMinObservationId()
        self.downloadTpfData()
        
        localUrl=self.dvFileLoc
        
        return self.parse(localUrl, *args, **kwargs)

