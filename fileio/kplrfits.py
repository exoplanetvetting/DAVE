# -*- coding: utf-8 -*-


#Note: binData is imported later in the file.
#This allows me a fall back to the slow code if it can't be found.
import numpy as np
#import nca
from metric_learn import nca


#Column definitions for the flux fits files. You no longer have to remember
#that PDC flux is in column 7, just ask for FluxColDef['PDCSAP_FLUX'].
#Not sure of the names? Try FLuxColDef.keys()
FluxColDef=dict(
             {   'TIME':0, \
                 'TIMECORR':1,
                 'CADENCENO':2, \
                 'SAP_FLUX':3, \
                 'SAP_FLUX_ERR':4,\
                 'SAP_BKG':5,\
                 'SAP_BKG_ERR':6,\
                 'PDCSAP_FLUX':7,\
                 'PDCSAP_FLUX_ERR':8,\
                 'SAP_QUALITY':9,\
                 'PSF_CENTR1':10,\
                 'PSF_CENTR1_ERR':11,\
                 'PSF_CENTR2':12,\
                 'PSF_CENTR2_ERR':13,\
                 'MOM_CENTR1':14,\
                 'MOM_CENTR1_ERR':15,\
                 'MOM_CENTR2':16,\
                 'MOM_CENTR2_ERR':17,\
                 'POS_CORR1':18,\
                 'POS_CORR2':19,\
             }
            )


SapQuality=dict(
            {
                'AttitudeTweak': 1, \
                'SafeMode':2, \
                'NotFinePoint':4, \
                'EarthPoint':8, \
                'ZeroCrossing':16, \
                'MomDump':32, \
                'BigArgabright':64, \
                'CosmicRay': 128, \
                'Exclude': 256, \
                'NotUsed': 512, \
                'SPSD': 1024, \
                'Outlier': 2048, \
                'SmallArgabright': 4096, \
                'CalibCosmic': 8192, \
                'LdeFlag': 16384, \
                'PdcCoarsePoint': 32768, \
                'NoData': 65536, \
                'RollingBandOptimal': np.uint32(2**17), \
                'RollingBandMask': np.uint32(2**18), \
                'PossibleRollTweak': np.uint32(2**19), \
                'DefiniteRollTweak': np.uint32(2**20), \
            }
          )


class KeplerQuarter():
    """A mapping from time (bkjd) <--> Quarter number

    Crude, but good enough
    """

    def __init__(self):
        quarterStart = np.zeros( (19,))
        quarterStart[0] = 120.538714567
        quarterStart[1] = 131.511972798
        quarterStart[2] = 169.764978159
        quarterStart[3] = 260.245273781
        quarterStart[4] = 352.397530014
        quarterStart[5] = 443.510594809
        quarterStart[6] = 539.470046009
        quarterStart[7] = 630.195605601
        quarterStart[8] = 735.384144927
        quarterStart[9] = 808.536246667
        quarterStart[10] = 906.866022230
        quarterStart[11] = 1001.228914371
        quarterStart[12] = 1099.429354807
        quarterStart[13] = 1182.757206809
        quarterStart[14] = 1274.159796224
        quarterStart[15] = 1373.508645059
        quarterStart[16] = 1472.117742730
        quarterStart[17] = 1559.246350576
        quarterStart[18] = 1591.001140092   #After end of last cadence

        #The numbers above are from a particular kepid (10214328)
        #I remove 15 mins to make sure the times are before BKJD for all
        #stars
        quarterStart[:18] -= 15*60/86400.

        #Add 1 cadences + slop to end time
        quarterStart[18] += 45*60/86400.

        self.quarterStart = quarterStart


    def getQuarterStartTime(self, quarter):
        if quarter < 0 or quarter > 17:
            raise ValueError("Invalid Quarter: %i" %(quarter))

        return self.quarterStart[quarter]

    def getQuarterForBkjd(self, bkjd):
        if bkjd < self.quarterStart[0]:
            raise ValueError("Requested time before start of Q0")

        if bkjd > self.quarterStart[18]:
            raise ValueError("Requested time after end of Q17")

        idx = self.quarterStart < bkjd
        return np.where(idx)[0][-1]


    def getQuartersForTransits(self, period_days, epoch_bkjd):
        """ Get an array of which quarter each transit falls in

        Inputs:
        ------
        period_days
            (float) Period of transit (days)

        epoch_bkjd
            (float) Time of first transit in BKJD


        Returns:
        --------
        An array. The ith element of the array is the quarter in which
        the ith transit falls.


        Notes:
        ---------
        If using this function to decide which quarters of data to load
        you should extract the unique elements of the list before
        passing the array to loadMultiQuarter()

        """
        per = float(period_days)
        t0 = float(epoch_bkjd)

        #Figure out when the transits are
        timespan = self.getQuarterTimeSpan(None)
        numTrans = np.ceil(timespan/per)
        i = np.arange(numTrans)
        trans = t0 + per*i
        trans = trans[trans < self.quarterStart[-1]]

        #Figure out which quarter each time is in.
        func = lambda x : np.where(x < self.quarterStart)[0][0]
        quarters = np.array(map(func, trans)) -1
        return quarters


    def getQuarterTimeSpan(self, quarter):
        if quarter is None:
            return self.quarterStart[-1] - self.quarterStart[0]

        else:
            return self.quarterStart[quarter+1] - self.quarterStart[quarter]


def loadMultiQuarterPdc(archive, kepid, quarterList, \
    fitCubic=True, removeEarthPoint=False):

    print("WARNING loadMultiQuarterPdc deprecated. Use")
    print("loadMultiQuarter() instead")
    return loadMultiQuarter(archive, kepid, quarterList, fitCubic=fitCubic, \
        removeEarthPoint=removeEarthPoint)


def loadMultiQuarter(archive, kepid, quarterList, \
    lcType='pdc', removeBadData=True, nCadencesForSmooth = None,
    fitCubic=False, removeEarthPoint=False, doMedianZero=True):
    """Loads many quarters of PDC and stitches them together

        Input:
        archive:        A mastio.KeplerArchive object
        kepid           (int) kepid of target of interest
        quarterList     (int list), which quarters to download

        Optional Inputs
        lcType                  Either PA or PDC. Case insensitive
        nCadencesForSmooth      If applying median filtering, how
                                wide should the filter be. Default is
                                to not smooth.
        fitCubic                Should a cubic polynomial be removed
                                from each quarter before merging.
                                Only useful for PDC.
        removeEarthPoint        Remove a number of cadences after
                                every cadence marked as Earthpoint.
        doMedianZero            Convert to fractional amplitude

        Returns:
        A two-d numpy array of all available data for the target for the
        requested quarter. Flagged cadences are removed. The PDC column
        is median zeroed and a cubic polynomial is removed from each
        quarter before stitching.

        Note:
        Fitting cubic is broken, and trying to use it will raise
        an exception (until I get around to fixing it)
    """

    if lcType.lower() =='pdc':
        fluxCol=FluxColDef['PDCSAP_FLUX']
    elif lcType.lower() =='pa':
        fluxCol=FluxColDef['SAP_FLUX']
    else:
        raise ValueError('lcType must be either PA or PDC')

    #import pdb; pdb.set_trace()
    #Remove any duplicate entries from quarterList
    quarterList = np.unique(quarterList)

    dataList = []
    for q in quarterList:
        #Load file
        try:
            fits = archive.getLongCadence(kepid, q)
        except IOError as e:
            print("WARN: Failed to load data for %i Q%i" %(int(kepid), q))
            print(e)
            continue

        data = getNumpyArrayFromFitsRec(fits)
        #Do this before we remove bad data (which includes
        #data with Earthpoint flag set)
        if removeEarthPoint:
            #import pdb; pdb.set_trace()
            idx = getIndicesOfEarthPointData(data)
            goodIdx = list(set(range(len(data))) - set(idx))
            data = data[goodIdx, :]

        #Remove bad data
        if removeBadData:
            idx = getIndicesOfGoodData(data)
            data = data[idx,:]

        #zero
        if doMedianZero:
            data = medianZero(data, fluxCol)

        if nCadencesForSmooth is not None:
            data = medianSubtract(data, nCadencesForSmooth, fluxCol)

        #Remove cubic trend
        if fitCubic:
            raise ValueError("Fitting cubic is broken. Don't try.")
            data = fitAndRemovePolynomial(data, 3)



        dataList.append(data)

    if len(dataList) == 0:
        raise IOError("Failed to find data in any listed quarter for %i" %(int(kepid)))

    out = np.concatenate( tuple(dataList))

    lookup = """ TIME TIMECORR CADENCENO
                 SAP_FLUX SAP_FLUX_ERR SAP_BKG SAP_BKG_ERR
                 PDCSAP_FLUX PDCSAP_FLUX_ERR SAP_QUALITY
                 PSF_CENTR1 PSF_CENTR1_ERR PSF_CENTR2 PSF_CENTR2_ERR
                 MOM_CENTR1 MOM_CENTR1_ERR MOM_CENTR2 MOM_CENTR2_ERR
                 POS_CORR1 POS_CORR2""".split()
    out = nca.Nca(out)
    out.setLookup(1, lookup)
    return out


def getNumpyArrayFromFitsRec(rec, cols=None):

    #If not defined, try to get all columns
    if cols is None:
        cols = range(len(rec.names))

    nRow = len(rec)
    nCol = len(cols)

    data = np.empty( (nRow, nCol) )

    for i in range(len(cols)):
        data[:,i] = rec.field(cols[i])

    return data


def getIndicesOfEarthPointData(data, numTrailingCadence=96):
    """Crude function to identify earthpoint cadences.

    It finds the earthpoint flag and flags the next 96
    cadences (2 days worth), regardless of their performance

    Inputs:
    data    (2d numpy array) PDC data

    Output:
    A list of indices of data possibly affected by earthpoints
    """


    flags = np.uint32(data[:, FluxColDef['SAP_QUALITY']])
    epFlag = np.uint32(SapQuality['EarthPoint'])
    ep = np.where(flags & epFlag)[0]


    idx = []
    for e in ep:
        tmp = np.arange(e, e+numTrailingCadence)
        idx.extend(tmp)

    return idx


def getIndicesOfBadData(data):
    allIdx = range(len(data))
    goodIdx = getIndicesOfGoodData(data)

    badIdx = set(allIdx) - set(goodIdx)
    return list(badIdx)


def getIndicesOfGoodData(data):
    """Remove data that is flagged as bad or has nans for time or flux

    Kepler data can be bad in three ways. The flux can be set to NaN,
    the *time* can be set to NaN, or a flag can be set. Not
    all flags are equally bad, so we include data flagged as
    zero crossing, spsd, or cosmic rays in the calibration data.
    These data may indicate lower quality data.


    Input:
    data:    (2d numpy array) Data to be smoothed, as loaded from loadData()

    Returns
    An array of ints containing the indices of good data
    """

    time = data[:, FluxColDef['TIME']]
    flux = data[:, FluxColDef['SAP_FLUX']]
    flag = data[:, FluxColDef['SAP_QUALITY']]

    #Set a mask that matches everything but listed flags.
    m = getMaskForBadData()
    mask = np.zeros( (len(flag)), dtype=np.uint32) + m

    #Find any cadences that match flag
    idx = np.bitwise_and(flag.astype(np.uint32), mask)
    idx = np.where(idx>0, 0, 1)    #Convert to booleans

    #Also match cadences where we don't flag but set either
    #time or flux to NaN
    idx = np.bitwise_and(idx, np.where(np.isfinite(time), 1, 0))
    idx = np.bitwise_and(idx, np.where(np.isfinite(flux), 1, 0))

    #Return only the indices set to true
    return np.where(idx)[0]


def getMaskForBadData():
    """If data[i, 'SAP_QUALITY] & m, data is marked with a flag
       indicating it's probably, or definiately bad.

       The flags listed here are considered informational only
       and not reasons to gap the cadence.
    """
    m = np.uint32(-1)
    m -= np.uint32(SapQuality['ZeroCrossing'])
    m -= np.uint32(SapQuality['SPSD'])
    m -= np.uint32(SapQuality['CosmicRay'])
    m -= np.uint32(SapQuality['CalibCosmic'])
    return m


def getMaskForBadK2Data():
    """If data[i, 'SAP_QUALITY] & m, data is marked with a flag
       indicating it's probably, or definiately bad.

    These are flags are appropriate for K2 data
    """
    m = np.uint32(0)
    m += np.uint32(SapQuality['AttitudeTweak'])
    m += np.uint32(SapQuality['NotFinePoint'])
    m += np.uint32(SapQuality['MomDump'])
    m += np.uint32(SapQuality['BigArgabright'])
    m += np.uint32(SapQuality['Exclude'])
    m += np.uint32(SapQuality['Outlier'])
    m += np.uint32(SapQuality['SmallArgabright'])
    m += np.uint32(SapQuality['NoData'])
    m += np.uint32(SapQuality['DefiniteRollTweak'])
    m += np.uint32(SapQuality['PossibleRollTweak'])

    return m


def medianSubtract(data, nPoints, fCol=FluxColDef['PDCSAP_FLUX']):
    """Remove trends from data

    Simplest possible lightcurve cleaning algorithm to clean data
    by subtracting off the median. Tested on single quarter data only,
    and enjoys only limited success.

    Inputs:
    data:    (2d numpy array) Data to be smoothed, as loaded from loadData()
    nPoints: How many adjacent points to be used in calculating the median
    fCol:    Which column of data contains the flux.

    Returns:
    data, with data[:,fCol] detrended
    """

    data[:,fCol] = medianSubtract1d(data[:, fCol], nPoints)
    return data


def medianSubtract1d(flux, nPoints):
    """Quick and dirty median smooth function.
    Not designed to be at all effecient"""

    size = len(flux)
    filtered = np.zeros_like(flux)
    for i in range(size):
        lwr = max(i-nPoints, 0)
        upr = min(lwr + 2*nPoints, size)
        lwr = upr- 2*nPoints

        sub = flux[lwr:upr]

        offset = np.median(sub)
        filtered[i] = flux[i] - offset

    return filtered


def fitAndRemovePolynomial(data, order, \
    timeCol=FluxColDef['TIME'], fluxCol=FluxColDef['PDCSAP_FLUX']):
    """

    Inputs:
    order   (int) order of fit. A good value for PDC of 3 (which fits
            cubic, i.e 4 parameters)
    """

    x = data[:,timeCol]
    y = data[:,fluxCol]
    polyCoeff = np.polyfit(x,y,order)
    bestFit = np.polyval(polyCoeff, x)

    data2 = data.copy() #Deep copy
    data2[:,fluxCol] -= bestFit
    return data2


def medianZero(data, col):
    """If x= data[:,col], set x = (x/median(x) - 1

        Note:
        This function is preferable to meanZero when the data has outliers.
    """

    med = np.median(data[:,col])

    tol=1e-9
    if np.fabs(med) < tol:
        raise ValueError("Median of lightcurve already nearly zero")

    data[:,col] /= med
    data[:,col] -= 1
    return data



###################################################








def old_foldAndBinData(x, y, period, epoch, expTime, numBins):

    assert(len(x) == len(y))
    dataSize = len(x)
    binWidth = period/float(numBins)

    value = np.zeros( (numBins) )
    weight = np.zeros( (numBins) )
    for i in range(dataSize):
        #Where is this cadence in binned space
        a0 = np.fmod(x[i] - epoch, period)
        a1 = np.fmod(x[i] - epoch + expTime, period)

        if a0 < 0:
            a0 += period
            a1 += period

        #Wrap a1 around if necessary
        if a1 < a0:
            a1 += period

        assert(0 <= a0 and a0 < period)
        assert(a0 <= a1)

        j = np.floor(a0/period * numBins)
        #Span (in bin space) of left most affected bin
        b0 = j*binWidth
        b1 = b0 + binWidth

        while b0 < a1:
            #How much does the cadence overlap with this bin
            overlap = min(a1, b1) - max(a0, b0)
            overlap /= a1-a0    #So overlap is a fraction E [0,1]

            value[j] += y[i]*overlap
            weight[j] += overlap

            #Move to the next bin and repeat
            b0 = b1
            b1 = b0+binWidth
            j = (j+1) % numBins

    #Normalise the values for the weights
    for i in range(numBins):
        value[i] /= weight[i]

    #Package the result for returning
    out = np.zeros( (numBins, 3) )
    #out[:,0] = np.arange(numBins)/float(numBins) * period + expTime/2.
    out[:,0] = np.arange(numBins)/float(numBins) * period
    out[:,1] = value
    out[:,2] = weight
    return out

try:
    from binData import foldAndBinData
except ImportError:
    print("INFO: binData.so not found. Falling back on old_foldAndBinnedData")
    foldAndBinData = old_foldAndBinData



def markTransitCadences(time, period_days, epoch_bkjd, duration_days,\
    numberOfDurations=1, flags=None):
    """Create a logical array indicating which cadences are
    affected by a transit

    Input:
    -------------
    time:
        (numpy 1d array) array of cadence times
    period_days:
        Transit period
    epoch_bkjd:
        Transit epoch
    duration_days:
        Duration of transit (start to finish). If you
        select a duration from first to last contact,
        all cadences affected by transit are selected.
        If you only select 2 to 3 contact, only the
        interior region of the transit is selected
    numberOfDurations
        How much of the lightcurve either side of
        the transit to mark. Default means to mark
        1/2 the transit duration either side of
        transit center.

    Optional Inputs:
    flags:
        If set, must be an array of booleans of length ``time``.
        Cadences where this array is true are ignored in the
        calculation. This is useful if some of the entries of time
        are Nan.

    Returns:
    -------------
    Array of length len(time) of booleans. Element set to true
    are affected by transits
    """

    if flags is None:
        flags = np.zeros_like(time, dtype=bool)

    i0 = np.floor((np.min(time[~flags])-epoch_bkjd)/period_days)
    i1 =  np.ceil((np.max(time[~flags])-epoch_bkjd)/period_days)
    assert(np.isfinite(i0))
    assert(np.isfinite(i1))

    irange = np.arange(i0, i1+1)
    transitTimes = epoch_bkjd + period_days*irange

    idx = np.zeros( len(time), dtype=np.bool8)
    for tt in transitTimes:
        diff = time - tt
        diff[flags] = 1e99  #A large value that isn't Nan
        assert(np.all(np.isfinite(diff)), "Nan found in diff")

        idx = np.bitwise_or(idx, np.fabs(diff) < \
            .5*duration_days*numberOfDurations)

    if not np.any(idx):
        print("WARN: No cadences found matching transit locations")
    return idx



def regularSpaceKeplerDataByCadence(data, missingValue=0,
    cCol=FluxColDef['CADENCENO']):
    """Take Kepler multi quarter data and insert zeros where no data is present.
    This results in regularly spaced data.

    Inputs:
    data    Numpy array as returned by kplrfits.loadMultiQuarter()

    Optional Inputs:
    missingValue    (float) What value to enter for cadences where there is
                    no data

    """

    nR,nCol = data.shape

    cin = data[:, cCol].astype(np.int32)
    nCadence = cin[-1]

    reg = np.zeros( (nCadence+1, nCol) ) + missingValue

    for i in range(nCol):
        reg[cin, i] = data[:,i]

    #Insert all the cadence values, even where data was missing
    reg[:, cCol] = np.arange(cin[0], cin[0]+nCadence+1)

    #Fix the infinities
    idx = np.isfinite(reg)
    reg[~idx] = missingValue

    reg = nca.Nca(reg, data.lookup)
    return reg
