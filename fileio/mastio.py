
__version__ = "$Id: mastio.py 1780 2014-08-27 16:36:11Z fmullall $"
__URL__ = "$URL: svn+ssh://flux/home/fmullall/svn/kepler/py/mastio.py $"


from dave.fileio.AbstractMast import KeplerAbstractClass
from dave.fileio.AbstractMast_TESS import TESSAbstractClass
from astroquery.mast import Observations
import numpy as np
import math
import os

"""
This module communicates with the MAST archive for the purpose of
downloading Kepler data. To reduce the amount of code written it uses
a fairly complex inheritance structure.

Summary
------------
To get K2 data do the following:
>>> ar = mastio.K2Archive()
>>> fits = ar.getLongCadence(epic, campaign)

To get the Vanderburg lightcurves, you use a VanderburgArchive class
>>> ar = mastio.VanderburgArchive()
>>> data = ar.getLongCadence(epic, campaign)


List of classes
-----------------
KeplerArchive
    Get data for classic data. Returns long and short cadence, lightcurves and TPF files

K2Archive
    Get official project data for K2. Returns all kinds of data. Lightcurves are
    processed with PDC

VanderburgArchive
    Lightcurves reduced with SFF. Long cadence lightcurves only. Data is returned as
    a numpy array

K2SC
    Lightcurves reduced by Suzanne Aigrain's method

Everest
    Lightcurves reduced by the Everest method.


Notes:
-----------
K2Archive() is used to get K2 data, KeplerArchive for classic Kepler data. Both
draw most of their code from AbstractKeplerClass() where all the common code is
written. Anything in the Abstract class can be used to query either archive.

To be even more general, functions to download and cache data from MAST are
implemented in the parent class of AbstractKeplerClass, called MastArchive. This
will make future efforts to get non-Kepler data from Mast a little easier.

Not every K2 reduction returns short cadence, or TPF files. All return
long cadence data.

Each class's getLongCadence() returns a different object depending on what
is stored at MAST.

Any optional arguments to getLongCadence() are passed directly to pyfits.getdata(),
so you can treat one function as a replacement for the other. K2Archive() also
has methods to get short cadence data, and target pixel files.
"""




##################################################################

class KeplerArchive(KeplerAbstractClass):
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

    def __init__(self, localDir=None, remoteServer="http://archive.stsci.edu"):
        remoteFluxPath = "pub/kepler/lightcurves"
        remoteTpfPath =  "pub/kepler/target_pixel_files"

        if localDir is None:
            localDir = os.path.join(os.environ['HOME'], ".mastio", "kepler")

        KeplerAbstractClass.__init__(self, localDir, remoteServer, remoteFluxPath, remoteTpfPath)

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


    def makeRemoteUrl(self, remotePath, kepid, quarter, filename, compressed):
        """Input arg **quarter** is not used in this function"""

        kidStr = "%09i" % (kepid)
        subdir = kidStr[:4]
        url = "%s/%s/%s/%s/%s" %(self.remoteServer, remotePath, subdir, \
                kidStr, filename)

        if compressed:
            url = url + '.gz'
        return url



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




#############################################################################
#############################################################################
#############################################################################


class K2Archive(KeplerAbstractClass):
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
        KeplerAbstractClass.__init__(self, localDir, remoteServer, remoteFluxPath, remoteTpfPath)


    def getFile(self, kepid, campaign, isFluxFile, month=None, *args, **kwargs):
        """See KeplerAbstractClass.getFile()"""

        #Traps the case of attempting to download a lightcurve file for the
        #first three quarters, for which no lightcurve files were created.
        if campaign < 1 and isFluxFile:
            raise ValueError("No lightcurve files were created for Campaigns 0,1 or 2")

        return KeplerAbstractClass.getFile(self, kepid, campaign, isFluxFile,
                                                month, *args, **kwargs)


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

#        if int(campaign) < 3:
#            subdir3 = "%i" %( 1e3*math.floor( subdir3/1e3))
#        else:
#            subdir3 = "%05i" %( 1e3*math.floor( subdir3/1e3))

        subdir3 = "%05i" % int(1e3*math.floor( subdir3/1e3))

        url = "%s/%s/%s/%s/%s/%s" %(self.remoteServer, remotePath, subdir1, subdir2, subdir3, filename)

        if compressed:
            url = url + '.gz'
        
        return url


#############################################################################
#############################################################################
#############################################################################

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
            raise ValueError("Vanderburg doesn't produce short cadence files")

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
        subdir3 = "%05i" %(int(kidStr[-5:]))

        #Note sure this is necessary
#        if int(campaign) < 3:
#            subdir3 = "%i" %( 1e3*math.floor( subdir3/1e3))

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


#############################################################################
#############################################################################
#############################################################################

class LocalK2Archive(K2Archive):
    def __init__(self, llcPath=".", lpdPath=None, slcPath=None, spdPath=None):
        """Get data from the local directories.

        A special case class to deal with a situation where data is available on
        a local disk in a directory structure different than what a K2Archive()
        wants to use. Each file type (long cadence, short cadence, flux/TPF file)
        is stored in its own directory.

        Optional Inputs:
        -----------
        llcPath
            (string) Location of long cadence flux files
        lpdPath
            (string) Location of long cadence TPF files. Defaults to llcPath
        slcPath
            (string) Location of short cadence flux files. Defaults to llcPath
        spdPath
            (string) Location of short  cadence TPF files. Defaults to **lpdPath**


        Notes:
        --------
        Does not cache files, nor does it seek to download the data if its missing.
       """

        #By setting remoteServer to a nonsensical value, we ensure we never
        #succeed in downloading any data
        K2Archive.__init__(self, remoteServer="www.example.com")
        self.llcPath = llcPath

        if lpdPath is None:
            lpdPath = self.llcPath
        self.lpdPath = lpdPath

        if slcPath is None:
            slcPath = self.llcPath
        self.slcPath = slcPath

        if spdPath is None:
            spdPath = self.lpdPath
        self.spdPath = spdPath


    def makeRemoteUrl(self, remotePath, kepid, quarter, filename, compressedOnServer):
        """This special case class should never go looking for
        remote files. If it does, it should definitely never find them"""
        return  "NoUrlForLocalArchiveClass"



#############################################################################
#############################################################################
#############################################################################

class K2SCArchive(KeplerAbstractClass):
    """Get co-trended long-cadence lightcurves as produced by Andrew Vanderburg
    from MAST.

    Vanderburg does not produce short cadence data, not does he
    make anything new available for TPF files, so this class can only
    return long cadence data.
    """

    def __init__(self, localDir = None, remoteServer="http://archive.stsci.edu"):
        if localDir is None:
            localDir = os.path.join(os.environ['HOME'], '.mastio', 'k2sc')

        remoteFluxPath = 'missions/hlsp/k2sc/v2'
        remoteTpfPath = 'NoK2ScTpfFiles'  #No TPFs for K2Sc
        KeplerAbstractClass.__init__(self, localDir, remoteServer, remoteFluxPath, remoteTpfPath)


    def getFilename(self, epic, campaign, isFluxFile, isShortCadence=False):
        """Vanderburg only creates LLC file equivalents.
        Rather than fix getFile(), to account for this, I know
        that getFile() always calls this method, so I can catch
        any attempts to get non-existent files here
        """
        if isShortCadence:
            raise ValueError("K2SC doesn't produce short cadence files")

        if not isFluxFile:
            raise ValueError("K2SC doesn't produce TPF files")

        version = 2
        #hlsp_k2sc_k2_llc_200004466-c03_kepler_v1_lc.fits
        fn = "hlsp_k2sc_k2_llc_%09i-c%02i_kepler_v%i_lc.fits"\
                %(epic, campaign, version)
        return fn


    def makeRemoteUrl(self, remotePath, k2id, campaign, filename, compressed):
        subdir1 = "c%02i" %(campaign)
        subdir2 = "%i" %(1e5*math.floor( int(k2id)/1e5))

        url = "%s/%s/%s/%s/%s" %(self.remoteServer, remotePath, \
                subdir1, subdir2, filename)


        if compressed:
            url = url + '.gz'
        return url



#############################################################################
#############################################################################
#############################################################################

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

        remoteFluxPath = 'missions/hlsp/everest/v2'
        remoteTpfPath = 'NoEverestTpfFiles'  #No TPFs for Everest
        KeplerAbstractClass.__init__(self, localDir, remoteServer, remoteFluxPath, remoteTpfPath)


    def getFilename(self, epic, campaign, isFluxFile, isShortCadence=False):
        """Vanderburg only creates LLC file equivalents.
        Rather than fix getFile(), to account for this, I know
        that getFile() always calls this method, so I can catch
        any attempts to get non-existent files here
        """
        if isShortCadence:
            raise ValueError("Everest doesn't produce short cadence files")

        if not isFluxFile:
            raise ValueError("Everest doesn't produce TPF files")

        version = 2 
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
        return url


#############################################################################
#############################################################################
#############################################################################

class TESSArchive(TESSAbstractClass):
    def __init__(self, localDir=None, remoteServer="http://archive.stsci.edu"):
#        """Download K2 target pixel files from the MAST archive
#
#        Input args:
#        localDir   (str) Where to cache data
#
#        Optional Input Args:
#        remoteServer:   (str) If MAST ever change the location of
#                        their server, the new url can be passed.
#
#        Usage:
#        #Create an archive object
#        archive = TESSArchive("/large_disk/TESS")
#
#        tpf, hdr = archive.getTargetPixelFile(k2id, c, header=True)
#        sc = archive.getTargetPixelFile(k2id, c, sc=True)
#
#
#        Note:
#        Only target pixel files are available from MAST at the moment
#        The interface is a little different from KeplerArchive()
#        """

#        if localDir is None:
#            localDir = os.path.join(os.environ['HOME'], '.mastio', 'tess')
        localDir = os.path.join(os.environ['HOME'], '.mastio', 'tess')

        remoteFluxPath = 'hlsps/tess-data-alerts'
        remoteTpfPath = 'hlsps/tess-data-alerts'
        TESSAbstractClass.__init__(self, localDir, remoteServer, remoteFluxPath, remoteTpfPath)


#https://archive.stsci.edu/hlsps/tess-data-alerts/hlsp_tess-data-alerts_tess_phot_00281459670-s01_tess_v1_tp.fits

    def getFile(self, kepid, campaign, isFluxFile, month=None, *args, **kwargs):
        """See KeplerAbstractClass.getFile()"""
        #Traps the case of attempting to download a lightcurve file for the
        #first three quarters, for which no lightcurve files were created.
        if campaign < 1 and isFluxFile:
            raise ValueError("No lightcurve files were created for Campaigns 0,1 or 2")

        return TESSAbstractClass.getFile(self, kepid, campaign, isFluxFile,
                                                month, *args, **kwargs)

    def getFilename(self, toi, campaign, isFluxFile, isShortCadence=False):

#	print isFluxFile

        if not isShortCadence:
            #Long cadence
            fType = 'tp'
            if isFluxFile:
                fType = 'lc'
        else:
            fType = 'spd-targ'
            if isFluxFile:
                fType = 'slc'

        version = 1
        fn = "hlsp_tess-data-alerts_tess_phot_%011i-s%02i_tess_v%1i_%s.fits" %(toi, campaign, version, fType)
        
        return fn

    def makeRemoteUrl(self, remotePath, toi, campaign, filename, compressed):

        url = "%s/%s/%s" %(self.remoteServer, remotePath, filename)

#        if compressed:
#            url = url + '.gz'
        return url
#############################################################################
#############################################################################
#############################################################################

class TESSArchive_AllData(TESSAbstractClass):
    def __init__(self, localDir=None, remoteServer="https://mast.stsci.edu"):#"https://archive.stsci.edu/"):#
#        """Download K2 target pixel files from the MAST archive
#
#        Input args:
#        localDir   (str) Where to cache data
#
#        Optional Input Args:
#        remoteServer:   (str) If MAST ever change the location of
#                        their server, the new url can be passed.
#
#        Usage:
#        #Create an archive object
#        archive = TESSArchive("/large_disk/TESS")
#
#        tpf, hdr = archive.getTargetPixelFile(k2id, c, header=True)
#        sc = archive.getTargetPixelFile(k2id, c, sc=True)
#
#
#        Note:
#        Only target pixel files are available from MAST at the moment
#        The interface is a little different from KeplerArchive()
#        """

#        if localDir is None:
#            localDir = os.path.join(os.environ['HOME'], '.mastio', 'tess')
        localDir = os.path.join(os.environ['HOME'], '.mastio', 'tess')

        remoteFluxPath = 'api/v0.1/Download/file/?uri=mast:TESS/product'#'missions/tess/tid/'#'hlsps/tess-data-alerts'
        remoteTpfPath = 'api/v0.1/Download/file/?uri=mast:TESS/product'#'missions/tess/tid/'#'hlsps/tess-data-alerts'

        TESSAbstractClass.__init__(self, localDir, remoteServer, remoteFluxPath, remoteTpfPath)

    def getFile(self, toi, sector, isFluxFile, month=None, *args, **kwargs):
        """See KeplerAbstractClass.getFile()"""
        #Traps the case of attempting to download a lightcurve file for the
        #first three quarters, for which no lightcurve files were created.
        if sector < 1 and isFluxFile:
            raise ValueError("No lightcurve files were created for Campaigns 0,1 or 2")

        return TESSAbstractClass.getFile(self, toi, sector, isFluxFile, month, *args, **kwargs)

    def getFilename(self, toi, sector, isFluxFile, isShortCadence=False):

        if not isShortCadence:
            #Long cadence
            fType = 'tp'
            if isFluxFile:
                fType = 'lc'
        else:
            fType = 'spd-targ'
            if isFluxFile:
                fType = 'slc'

        version = 1

#awkward TESS file naming
#
        toi_id_, sector_, tmp_ = "%016i" % toi, "s%04i/" % sector, "s%04i-%016i" % (sector, toi)
        prefix_ = "tess2018234235059-" if sector == 2 else "tess2018206045859-"
        suffix_ = "-0121-s_" if sector == 2 else "-0120-s_"

        fn = prefix_ + "s%04i-%016i" % (sector, toi) + suffix_+ "%s.fits" % fType 
#	fn_ = "s%04i/" % sector + toi_id_[0:4] + "/" + toi_id_[4:8] + "/" + toi_id_[8:12] + "/" + toi_id_[12:] + "/"
#	fn = fn_ + prefix_ + tmp_ + suffix_ + "%s.fits" % fType
#
#awkward TESS file naming
#	fn = "api/v0.1/Download/file/?uri=mast:TESS/product/tess*-s%04i-%016i*_%s.fits" %(sector, toi, fType)
#	fn = "s%04i/toi_id_[0:3]/toi_id_[4:7]/toi_id_[8:11]/toi_id_[12:-1]"
#	fn = "hlsp_tess-data-alerts_tess_phot_%011i-s%02i_tess_v%1i_%s.fits" %(toi, sector, version, fType)

#https://archive.stsci.edu/missions/tess/tid/s0002/0000/0000/3072/3470/tess2018234235059-s0002-0000000030723470-0121-s_lc.fits
#https://archive.stsci.edu/missions/tess/tid/s0002/0000/0003/0721/0830/tess2018234235059-s0002-0000000307210830-0121-s_tp.fits
#https://archive.stsci.edu/missions/tess/tid/s0001/0000/0003/0028/7949/tess2018206045859-s0001-0000000300287949-0120-s_lc.fits


        return fn

    def makeRemoteUrl(self, remotePath, toi, sector, filename, compressed):

        url = "%s/%s/%s" %(self.remoteServer, remotePath, filename)

#        if compressed:
#            url = url + '.gz'
        return url

