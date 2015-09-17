import urllib
import math
import os

import pyfits

__version__ = "$Id: mastio.py 1780 2014-08-27 16:36:11Z fmullall $"
__URL__ = "$URL: svn+ssh://flux/home/fmullall/svn/kepler/py/mastio.py $"


def createKeplerArchiveForLios():
    """Create an archive on my work machine, lios, caching data in the correct place"""
    localPath="/disk2/data"
    return KeplerArchive(localPath)

def createKeplerArchiveForUllord():
    localPath="/Users/fergal/data/kepler"
    return KeplerArchive(localPath)


def createKeplerArchiveForDunalmu():
    localPath="/home/fergal/data/fits/kepler"
    return KeplerArchive(localPath)


class MastArchive():
    """Base class for queries of the MAST Kepler and K2 archives.
    Handles submitting the URL, receiving and caching the results

    The child class determines the url to submit
    """
    def __init__(self, localDir, remoteServer, remoteFluxPath, remoteTpfPath):
        self.localDir = localDir
        self.remoteServer = remoteServer
        self.remoteFluxPath = remoteFluxPath
        self.remoteTpfPath = remoteTpfPath

        if not os.path.exists(localDir):
            raise ValueError("Local path %s doesn't exist" %(localDir))


    def getData(self, localUrl, remoteUrl, compressedOnServer, *args, **kwargs):

        forceDownload = kwargs.pop("force", False)

        if not forceDownload:
            if os.path.exists(localUrl):
                isEmpty = os.stat(localUrl)[6]==0
                if isEmpty:
                    raise IOError("File %s returns 404" %(remoteUrl))

        if not os.path.exists(localUrl):
            if compressedOnServer:
                self.download(remoteUrl, localUrl + '.gz')
                os.system("gunzip %s" %(localUrl + '.gz'))
            else:
                self.download(remoteUrl, localUrl)

        if 'ext' in kwargs and kwargs['ext'] == 0:
            try:
                retVal = pyfits.getheader(localUrl, *args, **kwargs)
            except IOError, e:
                raise e
        else:
            try:
                retVal = pyfits.getdata(localUrl, *args, **kwargs)
            except IOError, e:
                raise e

        return retVal


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
            return

        try:
            print remoteUrl
            (f,h) = urllib.urlretrieve(remoteUrl, localUrl)
#            print h
        except IOError, e:
            raise(e)


        if verbose:
            print h

        #try:
            #req = urllib2.Request(remoteUrl)
            ## create a request object

            #handle = urllib2.urlopen(req)
            ## and open it to return a handle on the url
        #except urllib2.HTTPError, e:
            #print 'We failed with error code - %s.' % e.code
            #if e.code == 404

        #fp = open(localUrl, "w")
        #fp.write(handle.read())
        #fp.close()

    def getLocalDir(self, quarter, kepid, quarterPrefix="Q"):
        quarterDir = "%s%s" %(quarterPrefix, quarter)

        localDir = os.path.join(self.localDir, quarterDir)
        if not os.path.exists(localDir):
            try:
                os.mkdir(localDir)
            except OSError, e:
                #Re-raising an exception makes the error easier to read
                raise e

        thou = "%04i" %( math.floor( float(kepid)/1000.))
        localDir = os.path.join(localDir, thou)
        if not os.path.exists(localDir):
            try:
                os.makedirs(localDir)
            except OSError. e:
                #Re-raising an exception makes the error easier to read
                raise e

        return localDir





class KeplerArchive(MastArchive):
    """Download public long cadence data from the classic Kepler
     archive. Keep a cached copy locally.

       Usage:
       1. Create a KeplerArchive object:
       localPath=/data/archive/kepler
       > archive = KeplerArchive(localPath)

       2. Retrieve a long cadence lightcurve on a target for a given
          quarter
       > kepid = 8112039
       > quarter = 5
       > data = archive.getLongCadence(kepid, quarter)

       The archive looks for a local copy of the lightcurve, and if
       one is not present, downloads it from the MAST website. It
       then loads the data using pyfits' getdata() command.

       All optional arguments to pyfits.getdata() are optional arguments
       to getLongCadence(). E.g.

       data, header = archive.getLongCadence(kepid, quarter, header=True)
       aperture = archive.getLongCadence(kepid, quarter, extname='aperture')


        Cached data is stored in localPath/Qn/ so data from different
        sessions can be reused.
    """

    def __init__(self, localDir=None, \
        remoteServer="http://archive.stsci.edu"):
        remoteFluxPath = "pub/kepler/lightcurves"
        remoteTpfPath =  "pub/kepler/target_pixel_files"

        if localDir is None:
            localDir = os.path.join(os.environ['HOME'], ".mastio", "kepler")

        MastArchive.__init__(self, localDir, remoteServer, remoteFluxPath, remoteTpfPath)

        self.timeStrs=dict([
            ('0', '2009131105131'), \
            ('1', '2009166043257'), \
            ('2', '2009259160929'), \
            ('3', '2009350155506'), \
            ('4', '2010078095331'), \
            ('5', '2010174085026'), \
            ('6', '2010265121752'), \
            ('7', '2010355172524'), \
            ('8', '2011073133259'), \
            ('9', '2011177032512'), \
            ('10', '2011271113734'), \
            ('11', '2012004120508'), \
            ('12', '2012088054726'), \
            ('13', '2012179063303'), \
            ('14', '2012277125453'), \
            ('15', '2013011073258'), \
            ('16', '2013098041711'), \
            ('17', '2013131215648'), \
        ])

        #Incomplete
        self.shortCadenceTimeStr=dict([
             ('2.1', '2009201121230'), \
             ('2.2', '2009231120729'), \
             ('5.1', '2010111051353'), \
             ('5.2', '2010140023957'), \
             ('5.3', '2010174090439'), \
             ('6.1', '2010203174610'), \
             ('6.2', '2010234115140'), \
             ('6.3', '2010265121752'), \
             ('7.1', '2010296114515'), \
             ('7.2', '2010326094124'), \
             ('7.3', '2010355172524'), \
             ('8.1', '2011024051157'), \
             ('8.2', '2011053090032'), \
             ('8.3', '2011073133259'), \
             ('9.1', '2011116030358'), \
             ('9.2', '2011145075126'), \
             ('9.3', '2011177032512'), \
            ('10.1', '2011208035123'), \
            ('10.2', '2011240104155'), \
            ('10.3', '2011271113734'), \
            ('11.1', '2011303113607'), \
            ('11.2', '2011334093404'), \
            ('11.3', '2012004120508'), \
            ('12.1', '2012032013838'), \
            ('12.2', '2012060035710'), \
            ('12.3', '2012088054726'), \
            ('13.1', '2012121044856'), \
            ('13.2', '2012151031540'), \
            ('13.3', '2012179063303'), \
            ('14.1', '2012211050319'), \
            ('14.2', '2012242122129'), \
            ('14.3', '2012277125453'), \
            ('15.1', '2012310112549'), \
            ('15.2', '2013011073258'), \
            ('15.3', '2012341132017'), \
            ('16.1', '2013017113907'), \
            ('16.2', '2013065031647'), \
            ('16.3', '2013098041711'), \
            ('17.1', '2013121191144'), \
            ('17.2', '2013131215648'), \
                   ] )

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
            compressedOnServer = True

        remoteUrl = self.makeRemoteUrl(remotePath, kepid, filename, compressedOnServer)
        return self.getData(localUrl, remoteUrl, compressedOnServer, *args, **kwargs)

##
    def getFilename(self, kepid, quarter, isFluxFile, month=None):
        if month is None:
            #Long cadence
            key = "%i" %(quarter)
            timestr = self.timeStrs[key]

            fType = 'lpd-targ'
            if isFluxFile:
                fType = 'llc'
        else:
            #Short cadence
            key = "%i.%i" %(quarter, month)
            timestr = self.shortCadenceTimeStr[key]

            fType = 'spd-targ'
            if isFluxFile:
                fType = 'slc'

        return "kplr%09i-%s_%s.fits" % \
            (int(kepid), timestr, fType)


    def makeRemoteUrl(self, remotePath, kepid, filename, compressed):

        kidStr = "%09i" % (kepid)
        subdir = kidStr[:4]
        url = "%s/%s/%s/%s/%s" %(self.remoteServer, remotePath, subdir, \
                kidStr, filename)

        if compressed:
            url = url + '.gz'
        return url




class K2Archive(MastArchive):
    def __init__(self, localDir=None, remoteServer="http://archive.stsci.edu"):
        """Download K2 target pixel files from the MAST archive

        Input args:
        localDir   (str) Where to cache data

        Optional Input Args:
        remoteServer:   (str) If MAST ever change the location of
                        their server, the new url can be passed.

        Usage:
        #Create an archive object
        archive = K2Archive("/large_disk/K2")

        tpf, hdr = archive.getTargetPixelFile(k2id, c, header=True)
        sc = archive.getTargetPixelFile(k2id, c, sc=True)


        Note:
        Only target pixel files are available from MAST at the moment
        The interface is a little different from KeplerArchive()
        """

        if localDir is None:
            localDir = os.path.join(os.environ['HOME'], '.mastio', 'k2')

        remoteFluxPath = 'pub/k2/lightcurves'
        remoteTpfPath = 'pub/k2/target_pixel_files'
        MastArchive.__init__(self, localDir, remoteServer, remoteFluxPath, remoteTpfPath)

    def getLongCadence(self, kepid, campaign, *args, **kwargs):
        """Get a long cadence lightcurve file.

        A thin wrapper around getFile()
        """
        return self.getFile(kepid, campaign, True, False, *args, **kwargs)


    def getShortCadence(self, kepid, campaign, *args, **kwargs):
        """Get a short cadence lightcurve file.

        A thin wrapper around getFile()
        """
        return self.getFile(kepid, campaign, True, True, *args, **kwargs)


    def getLongTpf(self, kepid, campaign, *args, **kwargs):
        """Get a long cadence TPF file.

        A thin wrapper around getFile()
        """
        return self.getFile(kepid, campaign, False, False, *args, **kwargs)


    def getShortTpf(self, kepid, campaign, *args, **kwargs):
        """Get a short cadence TPF file.

        A thin wrapper around getFile()
        """
        return self.getFile(kepid, campaign, False, True, *args, **kwargs)


    def getFile(self, kepid, campaign, isFluxFile, isShortCadence=False, *args, **kwargs):
        """Get a file of a given type for a given kepid/campaign[/month]

        Inputs:
        ---------
        kepid   (int)
            Kepid of target of interest
        campaign (int)
            campaign of interest
        isFluxFile (bool)
            If True, a lightcurve file is returned, otherwise a target
            pixel file is returned

        Optional Inputs:
        ----------------
        isShortCadence (bool)
            if True, get short cadence data

        All other optional inputs are passed to ``pyfits.getdata()``

        Returns:
        -------------
        A FITs object, and possibly more, depending on the optional
        arguments to ``pyfits.getdata``
        """

        filename = self.getFilename(kepid, campaign, isFluxFile, isShortCadence)
        localUrl = os.path.join( self.getLocalDir(campaign, kepid), filename)

        #Make remote url
        if isFluxFile:
            remotePath = self.remoteFluxPath
            compressedOnServer = False
        else:
            remotePath = self.remoteTpfPath
            compressedOnServer = True

        remoteUrl = self.makeRemoteUrl(remotePath, kepid, campaign, filename, compressedOnServer)
        return self.getData(localUrl, remoteUrl, compressedOnServer, *args, **kwargs)


    def getFilename(self, kepid, campaign, isFluxFile, isShortCadence=False):
        if not isShortCadence:
            #Long cadence
            fType = 'lpd-targ'
            if isFluxFile:
                fType = 'llc'
        else:
            fType = 'spd-targ'
            if isFluxFile:
                fType = 'slc'

        return "ktwo%09i-c%02i_%s.fits" % \
            (int(kepid), int(campaign), fType)


    def makeRemoteUrl(self, remotePath, k2id, campaign, filename, compressed):

        #/pub/k2/target_pixel_files/c0/200000000/00000

        kidStr = "%09i" % (k2id)
        subdir1 = "c%i" %(campaign)
        subdir2 = "%i" %(1e5*math.floor( int(k2id)/1e5))
        subdir3 = int(kidStr[-5:])


        if int(campaign) < 3:
            subdir3 = "%i" %( 1e3*math.floor( subdir3/1e3))
        else:
            subdir3 = "%05i" %( 1e3*math.floor( subdir3/1e3))

        url = "%s/%s/%s/%s/%s/%s" %(self.remoteServer, remotePath, \
                subdir1, subdir2, subdir3, filename)


        if compressed:
            url = url + '.gz'
        return url


#
