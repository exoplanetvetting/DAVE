import urllib
import math
import os

#from dave.fileio.mastio import KeplerAbstractClass
import dave.fileio.mastio
KeplerAbstractClass = dave.fileio.mastio.KeplerAbstractClass

import numpy as np

__version__ = "$Id: mastio.py 1780 2014-08-27 16:36:11Z fmullall $"
__URL__ = "$URL: svn+ssh://flux/home/fmullall/svn/kepler/py/mastio.py $"



class EverestArchive(KeplerAbstractClass):
    """Get co-trended long-cadence lightcurves as produced by Luger2016
    from MAST. This approach uses a pixel time series decorrelation technique.
    
    Everest does not produce short cadence data, not does it
    make anything new available for TPF files, so this class can only
    return long cadence data.
    """

    def __init__(self, localDir = None, remoteServer="http://archive.stsci.edu"):
        if localDir is None:
            localDir = os.path.join(os.environ['HOME'], '.mastio', 'everest')

        remoteFluxPath = 'missions/hlsp/everest'
        remoteTpfPath = 'NoEverestTpfFiles'  #No TPFs for Everest
        KeplerAbstractClass.__init__(self, localDir, remoteServer, remoteFluxPath, remoteTpfPath)
    

    def getFilename(self, epic, campaign, isFluxFile, isShortCadence=False):
        """Vanderburg only creates LLC file equivalents.
        Rather than fix getFile(), to account for this, I know
        that getFile() always calls this method, so I can catch
        any attempts to get non-existent files here
        """
        if isShortCadence:
            raise ValuerError("Everest doesn't produce short cadence files")

        if not isFluxFile:
            raise ValueError("Everest doesn't produce TPF files")
        
        version = 1
        #eg hlsp_everest_k2_llc_206103150-c03_kepler_v1.0_lc.fits  
        fn = "hlsp_everest_k2_llc_%09i-c%02i_kepler_v%3.1f_lc.fits"\
                %(epic, campaign, version)
        return fn


    def makeRemoteUrl(self, remotePath, k2id, campaign, filename, compressed):

        kidStr = "%09i" % (k2id)
        subdir1 = "c%02i" %(campaign)
        subdir2 = "%i" %(1e5*math.floor( int(k2id)/1e5))
        subdir3 = "%05i" % (int(kidStr[-5:]))

        url = "%s/%s/%s/%s/%s/%s" %(self.remoteServer, remotePath, \
                subdir1, subdir2, subdir3, filename)


        if compressed:
            url = url + '.gz'
#        print url
        return url

