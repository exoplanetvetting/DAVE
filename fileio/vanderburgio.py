import urllib
import math
import os

#from dave.fileio.mastio import KeplerAbstractClass
import dave.fileio.mastio
KeplerAbstractClass = dave.fileio.mastio.KeplerAbstractClass

import numpy as np

__version__ = "$Id: mastio.py 1780 2014-08-27 16:36:11Z fmullall $"
__URL__ = "$URL: svn+ssh://flux/home/fmullall/svn/kepler/py/mastio.py $"



class VanderburgArchive(KeplerAbstractClass):
    """Get co-trended long-cadence lightcurves as produced by Andrew Vanderburg
    from MAST. 
    
    Vanderburg does not produce short cadence data, not does he
    make anything new available for TPF files, so this class can only
    return long cadence data.
    """

    def __init__(self, localDir = None, remoteServer="http://archive.stsci.edu"):
        if localDir is None:
            localDir = os.path.join(os.environ['HOME'], '.mastio', 'vanderburg')

        remoteFluxPath = 'missions/hlsp/k2sff'
        remoteTpfPath = 'NoVanderburgTpfFiles'  #No TPFs for Vanderburg
        KeplerAbstractClass.__init__(self, localDir, remoteServer, remoteFluxPath, remoteTpfPath)
    

    def getFilename(self, epic, campaign, isFluxFile, isShortCadence=False):
        """Vanderburg only creates LLC file equivalents.
        Rather than fix getFile(), to account for this, I know
        that getFile() always calls this method, so I can catch
        any attempts to get non-existent files here
        """
        if isShortCadence:
            raise ValuerError("Vanderburg doesn't produce short cadence files")

        if not isFluxFile:
            raise ValueError("Vanderburg doesn't produce TPF files")
        
        version = 1
        fn = "hlsp_k2sff_k2_lightcurve_%09i-c%02i_kepler_v%i_llc-default-aper.txt"\
                %(epic, campaign, version)
        return fn


    def makeRemoteUrl(self, remotePath, k2id, campaign, filename, compressed):

        #/pub/k2/target_pixel_files/c0/200000000/00000

        kidStr = "%09i" % (k2id)
        subdir1 = "c%02i" %(campaign)
        subdir2 = "%i" %(1e5*math.floor( int(k2id)/1e5))
        subdir3 = int(kidStr[-5:])

        #Note sure this is necessary
        if int(campaign) < 3:
            subdir3 = "%i" %( 1e3*math.floor( subdir3/1e3))

        url = "%s/%s/%s/%s/%s/%s" %(self.remoteServer, remotePath, \
                subdir1, subdir2, subdir3, filename)


        if compressed:
            url = url + '.gz'
#        print url
        return url


    def parse(self, localUrl, *args, **kwargs):
        """Override MastArchive.parse()
        
        Mast data is usually stored as fits files, but Vanderburg lightcurves
        are plain text.
        """

        return np.loadtxt(localUrl, delimiter=",", skiprows=2, usecols=(0,1))
