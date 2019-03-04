#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  7 14:36:02 2018

@author: smullally

kwargs include ext=1 and header=True
"""

from astroquery.mast import Observations
from astroquery.mast import Tesscut
from astroquery.mast import Catalogs
from astropy.coordinates import SkyCoord
from astropy.io import fits


class TessAstroqueryArchive(object):
    
    def __init__(self, cachePath="./"):
        
        self.cachePath=cachePath

    def parse(self, localUrl, *args, **kwargs):
        """Load data using favourite IO class
        
        Assumes that localUrl exists on disk. Default file format
        for MAST is FITS, but sub-classes can override this method
        for other file types as necessary
        
        All optional arguments are passed into astropy.io.fits.getdata
        
        """
        if 'ext' in kwargs and kwargs['ext'] == 0:
            try:
                retVal = fits.getheader(localUrl, *args, **kwargs)
            except IOError as  e:
                raise e
        else:
            try:
                retVal = fits.getdata(localUrl, *args, **kwargs)
            except IOError as  e:
                raise e

        return retVal
    
    
    def getPostageStampObservationId(self, sector, ticid):
        """
        Using Astroquery find the two minute obsid.
        ticid = integer tic id
        sector = integer sector number
        """
        
        obs_query_string = "tess*-s%04i-%016i*" % (sector, ticid)
        obsTable = Observations.query_criteria(mission="TESS", obs_id=obs_query_string)
        
        #@todo check only one element returned
        #self.obsid = obsTable['obsid']
        
        return obsTable['obsid']
   
    
    def getOneFileFromMastObs(self, ticid, sector, fileType, *args, **kwargs):
        """
        Generic funtion to retrive a particular file type from a
        CAOM observation at the MAST.
        Private Function --  Only used by this class.
        """
        
        if fileType not in ['LC','TP','DVT']:
            raise ValueError("Requested file Type (%s) not in allowed list (LC,TP,DVT)."\
                             % fileType)
        
        obsid = self.getPostageStampObservationId(sector, ticid)
        
        manifest=Observations.download_products(obsid,\
                                            productSubGroupDescription=[fileType],\
                                            mrp_only=False,extension="fits",
                                            download_dir=self.cachePath)
        
        localUrl=manifest['Local Path'][0]
        
        return self.parse(localUrl, *args, **kwargs)
    
    
    def getLightcurve(self, ticid, sector, *args, **kwargs):
        """
        Get the light curve file for a given TIC Id and Sector.
        
        All optional arguments are passed into `astropy.io.fits.getdata`
        including header = True        
        
        """

        return self.getOneFileFromMastObs(ticid,sector,'LC', *args, **kwargs)
        
        
       
    
    def getDvt(self, ticid, sector, *args, **kwargs):
        """
        Get the Data Validation time series file for a given TIC Id and Sector.
        
        All optional arguments are passed into `astropy.io.fits.getdata`
        including header = True  
        ext = planetNumber
        
        """

        return self.getOneFileFromMastObs(ticid,sector,'DVT', *args, **kwargs)
        
    
    def getTPF(self, ticid, sector, *args, **kwargs):
        """
        Get the Target Pixel File for a given TIC Id and Sector.
        
        All optional arguments are passed into `astropy.io.fits.getdata`
        including header = True
        """
        
        return self.getOneFileFromMastObs(ticid,sector,'TP', *args, **kwargs)


    def getRaDec(self, ticid):
        """
        Use astroquery catalogs to get the RA and Dec
        """
        
        #@todo check ticid is an integer.
        
        starName="TIC %i" % int(ticid)
        cat = Catalogs.query_criteria(ID=starName,catalog = "Tic")
        coord = SkyCoord(cat[0]['ra'], cat[0]['dec'], unit = "deg")
        
        return coord
    
    def checkCoordInSector(self, coord, sector):
        """
        Determine if there is data on the Sector requested for the target.
        """

        sectorTable = Tesscut.get_sectors(coord)
        
        if sector in sectorTable['sector']:
            sectorExists=True
        else:
            sectorExists=False
            
        return sectorExists
    
    def getTessCut(self):
        """
        Perform tesscut.
        Not Done, come back.
        """
        
        manifest = Tesscut.get_cutouts(self.coord, size)
        
        
    
    def getFfiCutout(self,ticid, sector, size, *args, **kwargs):
        """
        Ask for astroquery Tesscut target pixel file.
        sizew needs to be in pixels.
        """
        
        coord = self.getRaDec(ticid)
        if not self.checkCoordInSector(coord,sector):
            raise ValueError("Requested TIC not in Requested Sector.")
                
        return self.getTessCut(coord, size)