
"""
Abstract classes for interfacing with the Mast Kepler Data archives.
 
Used by mastio.py
"""

import urllib
import math
import os

import pyfits


class MastArchive():
    """Base class for queries of the MAST Kepler and K2 archives.
    Handles submitting the URL, receiving and caching the results

    The child class determines the url to submit
    """
    def __init__(self, localDir, remoteServer):
        self.localDir = localDir
        self.remoteServer = remoteServer

        if not os.path.exists(localDir):
            try:
                os.mkdir(localDir)
            except IOError as e:
                msg = "Local path %s does not exist and can't be created: %s" \
                    %(localDir, e)
                raise IOError(e)


    def getData(self, localUrl, remoteUrl, compressedOnServer, *args, **kwargs):
        forceDownload = kwargs.pop("force", False)

        if not forceDownload:
            if os.path.exists(localUrl):
                isEmpty = os.stat(localUrl)[6]==0
                if isEmpty:
                    raise IOError("File %s previously returned 404" %(remoteUrl))

        if not os.path.exists(localUrl):
#	    print compressedOnServer	
            if compressedOnServer:
                self.download(remoteUrl, localUrl + '.gz')
                os.system("gunzip %s" %(localUrl + '.gz'))
            else:
                self.download(remoteUrl, localUrl)

        return self.parse(localUrl, *args, **kwargs)

    
    def download(self, remoteUrl, localUrl, clobber=True, verbose=False):
        """Download data from MAST archive for a given filename and a
        given quarter.

        Inputs:
        clobber:    If true, allow locally existing files to be overwritten
        verbose:    If true, print the file meta-information
        """

        if not clobber:
            if os.path.exists(localUrl):
                raise IOError("File %s already exists. Refusing to overwrite" %(localUrl))

        u = urllib.urlopen(remoteUrl)
        if u.getcode() == 404:
            #Write a zero length file to indicate no file is present
            if localUrl[-3:] == ".gz":
                localUrl = localUrl[:-3]

            fp = open(localUrl, 'w')
            fp.write("")
            fp.close()
            raise ValueError("File not found at MAST: %s" %(remoteUrl))
            return

        try:
            (f,h) = urllib.urlretrieve(remoteUrl, localUrl)
        except IOError as e:
            raise(e)


        if verbose:
            print(h)
    
    def parse(self, localUrl, *args, **kwargs):
        """Load data using favourite IO class
        
        Assumes that localUrl exists on disk. Default file format
        for MAST is FITS, but sub-classes can override this method
        for other file types as necessary
        """
        if 'ext' in kwargs and kwargs['ext'] == 0:
            try:
                retVal = pyfits.getheader(localUrl, *args, **kwargs)
            except IOError as e:
                raise e
        else:
            try:
                retVal = pyfits.getdata(localUrl, *args, **kwargs)
            except IOError as e:
                raise e

        return retVal

class TESSAbstractClass(MastArchive):
    """Base class for classes interacting with Kepler/K2/TESS data.
    
    This class contains methods for querying for Lightcurve and 
    TPF files in either short or long cadence.
    
    Example daughter classes are KeplerArchive and K2Archive.
    Daughter classes must, at a minimum, reimplement getFilename()
    and makeRemoteUrl(). You may want to write your own init, too.
    
    """
    
    def __init__(self, localDir, remoteServer, remoteFluxPath, remoteTpfPath):

        MastArchive.__init__(self, localDir, remoteServer)
        self.remoteFluxPath = remoteFluxPath
        self.remoteTpfPath = remoteTpfPath
        
    def getLocalDir(self, quarter, kepid, quarterPrefix="S"):

        quarterDir = "%s%s" %(quarterPrefix, quarter)

        localDir = os.path.join(self.localDir, quarterDir)

        if not os.path.exists(localDir):
            try:
                os.mkdir(localDir)
            except OSError as e:
                #Re-raising an exception makes the error easier to read
                raise e

        thou = str(int(kepid))#"%016i" % kepid#"%04i" %( math.floor( float(kepid)/1000.))

        localDir = os.path.join(localDir, thou)
        if not os.path.exists(localDir):
            try:
                os.makedirs(localDir)
            except OSError as e:
                #Re-raising an exception makes the error easier to read
                raise e

        return localDir


    def getLongCadence(self, kepid, quarter, *args, **kwargs):
        """Get a long cadence lightcurve file.

        A thin wrapper around getFile()
        """
        return self.getFile(kepid, quarter, True, None, *args, **kwargs)


    def getShortCadence(self, kepid, quarter, month, *args, **kwargs):
        """Get a short cadence lightcurve file.

        A thin wrapper around getFile()
        """
        return self.getFile(kepid, quarter, True, month, *args, **kwargs)


    def getLongTpf(self, kepid, quarter, *args, **kwargs):
        """Get a long cadence TPF file.

        A thin wrapper around getFile()
        """
        return self.getFile(kepid, quarter, False, None, *args, **kwargs)


    def getShortTpf(self, kepid, quarter, month, *args, **kwargs):
        """Get a short cadence TPF file.

        A thin wrapper around getFile()
        """
        return self.getFile(kepid, quarter, False, month, *args, **kwargs)


    def getFile(self, kepid, quarter, isFluxFile, month=None, *args, **kwargs):
        """Get a file of a given type for a given kepid/quarter[/month]

        Inputs:
        ---------
        kepid   (int)
            Kepid of target of interest
        quarter (int)
            Quarter of interest
        isFluxFile (bool)
            If True, a lightcurve file is returned, otherwise a target
            pixel file is returned

        Optional Inputs:
        ----------------
        month   (int)
            If None, a long cadence file is returned, if a value
            of 1,2 or 3 is supplied, a short cadence file is returned.

        All other optional inputs are passed to ``pyfits.getdata()``

        Returns:
        -------------
        A FITs object, and possibly more, depending on the optional
        arguments to ``pyfits.getdata``
        """

        filename = self.getFilename(kepid, quarter, isFluxFile, month)
        localUrl = os.path.join( self.getLocalDir(quarter, kepid), filename)

        #Make remote url
        if isFluxFile:
            remotePath = self.remoteFluxPath
            compressedOnServer = False
        else:
            remotePath = self.remoteTpfPath
            compressedOnServer = False#True

        remoteUrl = self.makeRemoteUrl(remotePath, kepid, quarter, filename, compressedOnServer)
        return self.getData(localUrl, remoteUrl, compressedOnServer, *args, **kwargs)


    def makeRemoteUrl(self, remotePath, kepid, quarter, filename, compressedOnServer):
        """Figure out the URL where the desired file is stored.
        
        Returns a string, typically similar to
        http://archive.stsci.edu/$remotePath/.../$filename
        """
        raise NotImplementedError("Not Implemented in Base Class")

    def getFilename(self, kepid, quarter, isFluxFile, isShortCadence):
        """Figure out the name of the file you want to download
        
        (But not its path)
        """
        raise NotImplementedError("Not Implemented in Base Class")

