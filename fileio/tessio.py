
"""
Created on Tue Nov 27 20:37:41 2018

@author: fergal
"""

from __future__ import print_function
from __future__ import division

from dave.fileio.AbstractMast import MastArchive
from astroquery.mast import Tesscut
from astropy.coordinates import SkyCoord
import os

class TessDvtLocalArchive(MastArchive):
    """
    Shim class to load DV Timeseries data from a local disk
        """    
    def __init__(self, path):
        MastArchive.__init__(self, "/dev/null", 'www.example.com')

        self.path = path

        self.filePrefix = dict()

#        self.filePrefix[1] = "tess2018206190142-s0001-s0001"        
#        self.filePrefix[2] = "tess2018235142541-s0002-s0002"
#	self.filePrefix[1] = "tess2019033200935-s0008-s0008"
#	self.filePrefix[1] = "hlsp_eleanor_tess_ffi_tic"
        
        self.fileSuffix = dict()
#        self.fileSuffix[1] = "00106_dvt.fits"
#        self.fileSuffix[2] = "00109_dvt.fits"
        self.fileSuffix[1] = "00182_dvt.fits"
#	self.fileSuffix[1] = "_s01_tess_v0.1.8_lc.fits"
        
    def getDvt(self, tic, sector, *args, **kwargs):

        """
        ext= to get different extensions
        header=True to return header
        """

        prefix= "tess2019033200935-s0008-s0008"#self.filePrefix[sector]
        suffix = "00182_dvt.fits"#self.fileSuffix[sector]

        localUrl = "%s-%016i-%s" %(prefix, int(tic), suffix)
#	localUrl = "%s%i%s" %(prefix, int(tic), suffix)
#	print(localUrl)

        localUrl = os.path.join(self.path, localUrl)

        if not os.path.exists(localUrl):
            raise IOError("File not found: %s" %(localUrl))
            
        return self.parse(localUrl, *args, **kwargs)


    def getTPF(self, tic, sector, *args, **kwargs):

        """
        ext= to get different extensions
        header=True to return header
        """

        prefix= "tess2019032160000-s0008"#self.filePrefix[sector]
        suffix = "0136-s_tp.fits"#self.fileSuffix[sector]

        localUrl = "%s-%016i-%s" %(prefix, int(tic), suffix)
#	localUrl = "%s%i%s" %(prefix, int(tic), suffix)
#	print(localUrl)

        localUrl = os.path.join(self.path, localUrl)

        if not os.path.exists(localUrl):
            raise IOError("File not found: %s" %(localUrl))
            
        return self.parse(localUrl, *args, **kwargs)


    def getTessCut(self):
        """
        Perform tesscut.
        Not Done, come back.
        """
        
        manifest = Tesscut.get_cutouts(self.coord, size)
