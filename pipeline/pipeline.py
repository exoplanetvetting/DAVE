# -*- coding: utf-8 -*-
"""
Created on Tue Oct  6 15:54:57 2015

@author: fergal

$Id$
$URL$
"""

__version__ = "$Id$"
__URL__ = "$URL$"




import numpy as np

import dave.pipeline.clipboard as clipboard
import dave.fileio.kplrfits as kplrfits

try:
    import dave.lpp.calcLPPoctave as lpp
except ImportError:
    print "Warn: LPP can't be imported"

import dave.fileio.mastio as mastio
import dave.fileio.tpf as tpf
import dave.fileio.nca as nca
import task
import os


def runOne(k2id, config):

    print "WARN: This function is deprecated. See main.py instead."
    taskList = config['taskList']

    clip = clipboard.Clipboard()
    clip['config'] = config
    clip['value'] = k2id
    #clip['dataStorePath'] = config['dataStorePath']

    #Check that all the tasks are properly defined
    for t in taskList:
        f = eval(t)

    #Now run them.
    for t in taskList:
        f = eval(t)
        clip = f(clip)

    if 'vet' in clip.keys():                                             # Jeff edits here - See if vetting has been done yet
        if 'reasonForFail' in clip.vet.keys():                           # Make sure a value for reason to fail exists else program may crash
            if 'ODD_EVEN_DIFF' in clip['vet.reasonForFail']:             # If it failed due to odd-even, then
                clip['bls.period'] = 2*clip['bls.period']                # Double the period
                taskList = """trapezoidFitTask vetTask plotTask""".split()  # And re-fit, re-vet, and re-plot
                clip['vet.fluxVet.comments'] = clip['vet.fluxVet.comments'] + "Re-fit at twice period due to odd/even"   # Make a note we re-fit at 2*period
                for t in taskList:
                    f = eval(t)
                    clip = f(clip)

    print "WARN: This function is deprecated. See main.py instead."
    return clip


#
def loadDefaultConfig():
    cfg = dict()
    cfg['debug'] = True
    cfg['campaign'] = 3
    cfg['timeout_sec'] = 120
    cfg['nPointsForMedianSmooth'] = 48
    cfg['blsMinPeriod'] = 0.5
    cfg['blsMaxPeriod'] = 30


    #Thermal resettling tends to bugger up the first day or so of K2 data
    #It seems best to just cut out a chunk.
    cfg['numInitialCadencesToIgnore'] = 50

    #Vetting parameters
    cfg['minSnrForDetection'] = 10.
    cfg['maxLppForTransit'] = 10**-2.1
    #How significant must centroid offset be to claim a false positive.
    #This value is between 0 and 1. Higher values mean fewer false positives
    cfg['minProbDiffImgCentroidForFail'] = 0.99


    #The LPP mapping file is used by LPP to define the regions of param
    #space where the transits cluster.
    try:
        path = lpp.getLppDir()
        cfg['lppMapFilePath'] = os.path.join(path, "octave/maps/mapQ1Q17DR24-DVMed6084.mat")
    except NameError:
        cfg['lppMapFilePath'] = "NotNeeded"

    davePath = os.path.join(os.environ['HOME'],"daveOutput","")
    cfg['modshiftBasename'] =  davePath
    cfg['onepageBasename']  = davePath
    cfg['clipSavePath'] = davePath
    #Location of the place all the light curves and TPF files are stored
    cfg['dataStorePath'] = os.path.join(os.environ['HOME'],".mastio/k2")

    #Location of the model PRF fits files.
    cfg['prfPath'] = os.path.join(os.environ['HOME'], ".mastio/keplerprf")

    #My front end
    tasks = """checkDirExistTask serveTask extractLightcurveTask
        computeCentroidsTask rollPhaseTask cotrendSffDataTask detrendDataTask'
        fblsTask trapezoidFitTask vetTask plotTask saveClip""".split()   # Jeff added plotTask
    cfg['taskList'] = tasks

    searchTaskList = """blsTask trapezoidFitTask modshiftTask
    measureDiffImgCentroidsTask dispositionTask""".split()
    cfg['searchTaskList'] = searchTaskList

    cfg['keysToIgnoreWhenSaving'] = ["serve"]
    return cfg




@task.task
def checkDirExistTask(clip):
    """
    Code checks to see if certain directories or files exist.
    """
    modbase=clip['config.modshiftBasename']
    plotbase=clip['config.onepageBasename']  # Jeff added this
    prfdir=clip['config.prfPath']
    lppmap=clip['config.lppMapFilePath']
    datadir=clip['config.dataStorePath']
    clipdir=clip['config.clipSavePath']

    moddir=os.path.dirname(modbase)
    plotdir=os.path.dirname(plotbase)
    vetdir=getVetDir()
    vetExName="%s%s" % (vetdir, "/modshift")

    errors=[]
    try:
        open(lppmap)
    except IOError:
        errors.append(("Cannot Find LPP Map, %s " % lppmap))
    if not (os.path.exists(moddir)):
        errors.append("Cannot Find Modshift Write Dir, %s " % moddir)
    if not (os.access(moddir, os.W_OK)):
        errors.append("Cannot Write to modshift directory %s " % moddir)
    if not (os.path.exists(plotdir)):
        errors.append("Cannot Find Plotting Write Dir, %s " % plotdir)
    if not (os.access(plotdir, os.W_OK)):
        errors.append("Cannot Write to plotting directory %s " % plotdir)
    if not (os.access(clipdir,os.W_OK)):
        errors.append("Can not write to clip directory %s" % clipdir )
    if not (os.access(prfdir,os.R_OK)):
        errors.append("Cannot Read from prf Directory, %s "% prfdir)
    if not (os.access(datadir,os.R_OK)):
        errors.append("Cannot Read from Data Dir %s " % datadir)
    if not (os.access(vetExName,os.EX_OK)):
        errors.append("%s not found. Did you run make?" % vetExName )


    if len(errors) > 0:
        new="\n"
        msg=new.join(errors)
        raise IOError(msg)

    return clip


def getVetDir():
    """Get the path where Vetting stores its executables"""
    pathSep = "/"
    path = os.path.realpath(__file__)
    path = pathSep.join(path.split(pathSep)[:-2])
    path = "%s%s" % (path,"/vetting")
    return path


@task.task
def serveTask(clip):
    k2id = clip['value']
    campaign = clip['config.campaign']
    storeDir = clip['config.dataStorePath']
    detrendType = clip.get('config.detrendType', 'PDC')

    ar = mastio.K2Archive(storeDir)
    clip['serve'] = loadTpfAndLc(k2id, campaign, ar, detrendType)

    #Enforce contract. (Make sure expected keys are in place)
    clip['serve.time']
    clip['serve.cube']
    clip['serve.socData']
    clip['serve.tpfHeader']
    return clip

@task.task
def serveLocalTask(clip):
    k2id= clip['value']
    campaign = clip['config.campaign']
#    lcpath="/soc/nfs/production-nfs2/c7/exports/archive_ksop2554/lcv/"
#    tppath="/soc/nfs/production-nfs2/c7/exports/archive_ksop2554/cad_targ_soc/"
    lcpath="/media/soc-nfs/production_nfs2/c7/exports/archive_ksop2554/lcv/"
    tppath="/media/soc-nfs/production_nfs2/c7/exports/archive_ksop2554/cad_targ_soc/"

    ar = mastio.LocalK2Archive(llcPath=lcpath,lpdPath=tppath)
    clip['serve'] = loadTpfAndLc(k2id,campaign,ar)

    clip['serve.time']
    clip['serve.cube']
    clip['serve.socData']
    clip['serve.tpfHeader']
    return clip


@task.task
def extractLightcurveTask(clip):
    time = clip['serve.time']
    data = clip['serve.socData']
    numInitialCadencesToIgnore = clip['config.numInitialCadencesToIgnore']
    flagValues = clip.get('serve.flags', data[:, 'SAP_QUALITY'])
    flux = data[:, 'SAP_FLUX'].copy()

    #Convert flags to a boolean, and flag other bad data
    mask = kplrfits.getMaskForBadK2Data()
    flags = (flagValues & mask).astype(bool)
    flags |= ~np.isfinite(time)
    flags |= ~np.isfinite(flux)
    #flags[flux<1] = True
    flags[:numInitialCadencesToIgnore] = True

    #Placeholder. Use the SOC PA data for the lightcurve
    out = dict()
    out['rawLightcurve'] = flux
    out['time'] = time
    clip['extract'] = out
    clip['extract.source'] = "SOC PA Pipeline"
    clip['extract.flags'] = flags

    #Enforce contract
    clip['extract.rawLightcurve']
    clip['extract.flags']
    return clip


@task.task
def cotrendDataTask(clip):
    """Produce a cotrended lightcurve in units of fractional amplitude"""

    data = clip['serve.socData']
    flags = clip['extract.flags']
    flux = data[:, 'PDCSAP_FLUX'].copy()

    #Cotrending may also produce Nans
    flags |= ~np.isfinite(flux)

    #Remove dc offset
    dcOffset = np.median( flux[~flags])
    flux = (flux/ dcOffset) - 1
    clip['cotrend'] = {'flux_frac': flux}
    clip['cotrend.dcOffset'] = dcOffset
    clip['cotrend.flags'] = flags
    clip['cotrend.dcOffset'] = dcOffset
    clip['cotrend.source'] = "SOC PDC Pipeline"

    #Enforce contract
    clip['cotrend.flux_frac']
    return clip



from dave.blsCode import outlier_detection
import dave.misc.noise as noise
@task.task
def detrendDataTask(clip):
    """

    TODO:
    This code could be considerably simplified if we set flux[singleOutlierIndices]
    ==0. Would this produce the same results?
    """


    flux = clip['cotrend.flux_frac']
    time = clip['serve.time']
    flags = clip['cotrend.flags']

    nPoints = clip['config.nPointsForMedianSmooth']

    #When you detrend, you must do something about the gaps and bad values.
    #This is the simplest possible thing. Replace all bad/missing data with
    #zeros. This is a placehold. Bad data inside a transit is replaced with
    #a zero.
    #flux[flags] = 0
    #Flag the outliers in the data
    #singleOutlierIndices = outlier_detection.outlierRemoval(time, flux)
    flags |= noise.singlePointDifferenceSigmaClip(flux, initialClip=flags)

    #do a simple detrend with the outliers not included
    trend = outlier_detection.medianDetrend(flux[~flags], nPoints)
    medianDetrend = np.zeros_like(flux)
    medianDetrend[~flags] = flux[~flags] - trend
    

    flags |= (medianDetrend < -1)    

    clip['detrend'] = dict()
    clip['detrend.flux_frac'] = medianDetrend
    clip['detrend.flags'] = flags
    clip['detrend.source'] = "Simple Median detrend"
    
    return clip



@task.task
def computeCentroidsTask(clip):
    data = clip['serve.socData']

    cent_colrow = np.empty( (len(data), 2))
    cent_colrow[:,0] = data[:, 'MOM_CENTR1'].copy()
    cent_colrow[:,1] = data[:, 'MOM_CENTR2'].copy()
    clip['centroids'] = {'cent_colrow': cent_colrow}
    clip['centroids.source'] = "SOC PA Pipeline"

    #Enforce contract
    clip['centroids.cent_colrow']
    return clip


import dave.diffimg.arclen as arclen
@task.task
def rollPhaseTask(clip):

    centColRow = clip['centroids.cent_colrow']
    flags = clip['extract.flags']
    rot = arclen.computeArcLength(centColRow, flags>0)
    rollPhase = rot[:,0]
    rollPhase[flags>0] = -9999    #A bad value

    clip['rollPhase'] = {'rollPhase':rollPhase}
    return clip





@task.task
def singleEventSearchTask(clip):


    clip['eventList'] = []
    subClip = searchForEvent(clip)

    if 'exception' in subClip.keys():
        clip['exception'] = subClip['exception']
        clip['backtrace'] = subClip['backtrace']

    clip['eventList'].append(subClip)
    return clip


def searchForEvent(clip):
    subClip = clip.shallowCopy()

    originalKeyList = subClip.keys()
    taskList = clip['config.searchTaskList']

    #Set the flags attribute of the new subclip
    #Problem with this code is it closely tied to the behaviour
    #of multiEventSearchTask
    try:
        tmp = clip.eventList[-1]
        flags = tmp['flags']
    except (IndexError, KeyError):
        flags = clip['detrend.flags']
    subClip['flags'] = flags

    #Check that all the tasks are properly defined
    for t in taskList:
        f = eval(t)

    #Now run them.
    for t in taskList:
        f = eval(t)
        subClip = f(subClip)

#    #@TODO List of tasks to run should be config param
#    subClip = placeholderBls(subClip)
#    subClip = trapezoidFitTask(subClip)
#    subClip = modshiftTask(subClip)
#    subClip = measureDiffImgCentroidsTask(subClip)
#    subClip = dispositionTask(subClip)

    newKeys = list(set(subClip.keys()) - set(originalKeyList))
    out = clipboard.Clipboard(__meta__=subClip['__meta__'])
    for k in newKeys:
        out[k] = subClip[k]


    #Mark all locations for this event as data not to be used.
    time = subClip['serve.time']
    period_days = subClip['trapFit.period_days']
    epoch_bkjd = subClip['trapFit.epoch_bkjd']
    duration_days = subClip['trapFit.duration_hrs'] / 24.

#    assert(np.all(np.isfinite(time[~flags])))
#    assert(np.any(flags))
    idx = kplrfits.markTransitCadences(time, period_days, epoch_bkjd, \
        duration_days, numberOfDurations=2, flags=flags)

    out['flags'] = flags | idx

    return out


import dave.blsCode.fbls as fbls
import dave.misc.noise as noise
@task.task
def fblsTask(clip):
    time_days = clip['extract.time']
    flux_norm = clip['detrend.flux_frac']
    flags = clip['detrend.flags']
    minPeriod = clip['config.blsMinPeriod']
    maxPeriod = clip['config.blsMaxPeriod']

#    durations = np.array([ 2,4,6,8, 10, 12])/24.
    durations = np.array([ 4,6,8, 10, 12])/24.
    idx = flags == 0
    blsObj = fbls.BlsSearch(time_days[idx], flux_norm[idx], \
        [minPeriod, maxPeriod], durations)

    period, epoch, depth, duration = blsObj.getEvent()
    spectrum = blsObj.compute1dBls()

    duration_cadences = int(np.round(duration*48)) #Correct for K2
    rms = noise.computeSgCdpp_ppm(flux_norm[idx], duration_cadences)*1e-6
    idx = kplrfits.markTransitCadences(time_days, period, epoch, \
        duration, flags=flags)
    snr = (depth/rms)*np.sqrt(np.sum(idx))

    out = dict()
    out['period'] = period
    out['epoch'] = epoch
    out['duration_hrs'] = duration * 24
    out['depth'] = depth
    out['snr'] = snr
    out['bls_search_periods'] = spectrum[:,0]
    out['convolved_bls'] = spectrum[:,1]
    #out['obj'] = blsObj
#    out['bls'] = bls  #bls array is extremely big
    clip['bls'] = out

    #Enforce contract
    clip['bls.period']
    clip['bls.epoch']
    clip['bls.duration_hrs']
    return clip


#import dave.blsCode.bls_ktwo as bls
#@task.task
#def blsTask(clip):
    #time_days = clip['serve.time']
    #flux_norm = clip['detrend.flux_frac']
    #flags = clip['detrend.flags']
    #minPeriod = clip['config.blsMinPeriod']
    #maxPeriod = clip['config.blsMaxPeriod']

    ##Zero out the bad data. This crashes BLS
##    flux_norm[flags] = 0
##    assert(np.all( np.isfinite(flux_norm)))

    #idx = flags == 0
    #period, epoch, duration, depth, bls_search_periods, convolved_bls = \
        #bls.doSearch(time_days[idx], flux_norm[idx], minPeriod, maxPeriod)

    #out = clipboard.Clipboard()
    #out['period'] = period
    #out['epoch'] = epoch
    #out['duration_hrs'] = duration * 24
    #out['depth'] = depth
    #out['bls_search_periods'] = bls_search_periods
    #out['convolved_bls'] = convolved_bls
    #clip['bls'] = out

    ##Enforce contract
    #clip['bls.period']
    #clip['bls.epoch']
    #clip['bls.duration_hrs']
    #return clip


@task.task
def placeholderBls(clip):
    """Debugging code. Returns the ephemeris of the largest event in
    K2Id 206103150
    """
    out = clipboard.Clipboard()
    out['period'] = 4.15892
    out['epoch'] = 2145.76
    out['duration_hrs'] = 1.94443
    out['depth'] = .01112825

    clip['bls'] = out
    return clip


@task.task
def lppMetricTask(clip):
    time_days = clip['serve.time']
    flux_unitmean = clip['detrend.flux_frac']+1
    period_days = clip['bls.period']
    duration_hrs = clip['bls.duration_hrs']
    phase_bkjd = clip['bls.epoch']  #Check this what BLS returns
    mapFile = clip['config.lppMapFilePath']

    TLpp, Y, binnedFlux = lpp.fergalVersion(time_days, flux_unitmean, mapFile,\
        period_days, duration_hrs, phase_bkjd)

    out = dict()
    out['TLpp'] = TLpp
    out['binnedFlux'] = binnedFlux

    clip['lpp'] = out

    #Enforce contract
    clip['lpp.TLpp']

    return clip


import dave.trapezoidFit.estimateSnr as tf
@task.task
def trapezoidFitTask(clip):
    time_days = clip['serve.time']
    flux_norm = clip['detrend.flux_frac']
    flags = clip['detrend.flags']
    period_days = clip['bls.period']
    duration_hrs = clip['bls.duration_hrs']
    phase_bkjd = clip['bls.epoch']  #Check this what BLS returns
    depth_frac = np.abs(clip['bls.depth'])

    #We don't know these values.
    unc = np.ones_like(flux_norm)
    unc[flags] = 1e99
    flux_norm[flags] = 0


    assert(np.all(np.isfinite(time_days[~flags])))
    assert(np.all(np.isfinite(flux_norm[~flags])))
    out = tf.getSnrOfTransit(time_days, flux_norm,\
        unc, flags, \
        period_days, phase_bkjd, duration_hrs, depth_frac)

    assert(len(time_days) == len(out['bestFitModel']))
    clip['trapFit'] = out

    clip['trapFit.period_days']
    clip['trapFit.epoch_bkjd']
    clip['trapFit.duration_hrs']
    clip['trapFit.ingress_hrs']
    clip['trapFit.depth_frac']
    clip['trapFit.bestFitModel']
    clip['trapFit.snr']
    return clip



#import dave.trapezoidFit.trapfit as trapFit
import dave.vetting.ModShift as ModShift
@task.task
def modshiftTask(clip):

    time = clip['serve.time']
    flux = clip['detrend.flux_frac']
    fl = clip['detrend.flags']
    period_days = clip['trapFit.period_days']
    epoch_bkjd = clip['trapFit.epoch_bkjd']
    #dur_hrs =  clip['trapFit.duration_hrs']
    #ingress_hrs = clip['trapFit.ingress_hrs']
    #depth_ppm = 1e6*clip['trapFit.depth_frac']

    epic = clip['value']
    basePath = clip['config.modshiftBasename']
    epicStr = "%09i" %(epic)
    basename = getOutputBasename(basePath, epicStr)

    # Name that will go in title of modshift plot
    objectname = "EPIC %09i" %(epic)
#
#    subSampleN= 15
#    ioBlock = trapFit.trapezoid_model_onemodel(time[~fl], period_days, \
#        epoch_bkjd, depth_ppm, dur_hrs, \
#        ingress_hrs, subSampleN)
#    model = ioBlock.modellc -1   #Want mean of zero
#    #model *= -1  #Invert for testing
    model = clip['trapFit.bestFitModel']
    modplotint=1  # Change to 0 or anything besides 1 to not have modshift produce plot
    plotname = "%s-%02i-%04i-%s" % (basename,np.round(clip.bls.period*10),np.round(clip.bls.epoch),clip.config.detrendType)
    out = ModShift.runModShift(time[~fl], flux[~fl], model[~fl], \
        plotname, objectname, period_days, epoch_bkjd, modplotint)
    clip['modshift'] = out

    #Enforce contract
    clip['modshift.mod_Fred']
    clip['modshift.mod_ph_pri']
    clip['modshift.mod_secdepth']
    clip['modshift.mod_sig_pri']
    return clip


import dave.diffimg.centroid as cent
import dave.diffimg.prf as prf
@task.task
def measureDiffImgCentroidsTask(clip):

    #Measuring centroids requires a lot of input params
    period_days = clip['trapFit.period_days']
    epoch_bkjd = clip['trapFit.epoch_bkjd']  #Check this what BLS returns
    duration_hrs = clip['trapFit.duration_hrs']

    cube = clip['serve.cube']
    cube[ ~np.isfinite(cube) ] = 0
    tpfHeader0 = clip['serve.tpfHeader0']
    tpfHeader = clip['serve.tpfHeader']
    ccdMod = tpfHeader0['MODULE']
    ccdOut = tpfHeader0['OUTPUT']
    bbox = cent.getBoundingBoxForImage(cube[0], tpfHeader)
    rollPhase = clip['rollPhase.rollPhase']
    prfPath = clip['config.prfPath']
    prfObj = prf.KeplerPrf(prfPath)

    time_days = clip['serve.time']
    qflags = clip['serve.flags']
    flags = clip['detrend.flags']

#    import pdb; pdb.set_trace()
    out, diagnostics, log = cent.measureDiffOffset(period_days, epoch_bkjd, duration_hrs, \
        time_days, prfObj, ccdMod, ccdOut, cube, bbox, rollPhase, flags, qflags)

    #Set column names
    out = nca.Nca(out)
    out.setLookup(1, "rin intr_col intr_row diff_col diff_row".split())

    clip['diffImg'] = {'centroid_timeseries':out, \
                       'diagnostics':diagnostics, \
                       'log':log}

    clip['diffImg.centroid_timeseries']
    return clip


#import dave.vetting.RoboVet as RoboVet
@task.task
def vetTask(clip):

    print "WARN: vetTask is deprecated. Use dispostionTask instead"

    snr = clip['trapFit.snr']
    snrThreshold = clip['config.minSnrForDetection']
    lppThreshold = clip['config.maxLppForTransit']
    offsetThreshold_sigma = clip['config.minCentroidSignifForFp']

    out = clipboard.Clipboard(isCandidate=True, reasonForFail="None")
    if snr < snrThreshold:
        out['isCandidate'] = False
        out['reasonForFail'] = "SNR (%.1f) below threshold %.1f" \
            %(snr, snrThreshold)
#        clip['vet'] = out
#        return clip

    clip = lppMetricTask(clip)
    if 'exception' in clip.keys():
        return clip

#    Tlpp = clip['lpp.TLpp']
#    if Tlpp > lppThreshold:
#        out['isCandidate'] = False
#        out['reasonForFail'] = "TLpp (%.1f) above threshold %.1f" \
#            %(Tlpp, lppThreshold)
##        clip['vet'] = out
##        return clip

    clip = modshiftTask(clip)
    if 'exception' in clip.keys():
        return clip
    modshift = clip['modshift']
    fluxVetDict = RoboVet.roboVet(modshift)
    out['fluxVet'] = fluxVetDict

    assert(fluxVetDict['disp'] in ["candidate", "false positive"])

    if fluxVetDict['disp'] == "false positive":
        out['isCandidate'] = False
        out['reasonForFail'] = fluxVetDict['comments']
#        clip['vet'] = out
#        return clip

    clip = measureDiffImgCentroidsTask(clip)
    if 'exception' in clip.keys():
        return clip
    centroids = clip['diffImg.centroid_timeseries']

    result = cent.measureOffsetInTimeseries(centroids)
    out['centroidVet'] = result
    signif = result['signif']
    offset = result['offset']

    if signif > offsetThreshold_sigma:
        out['isCandidate'] = False
        out['reasonForFail'] = "Centroid offset of %.2f (%.1f sigma) detected" \
            %( offset, signif)
#        clip['vet'] = out
#        return clip


    clip['vet'] = out

    #Enforce contract
    clip['vet.isCandidate']
    clip['vet.reasonForFail']
    return clip



import dave.vetting.RoboVet as RoboVet
@task.task
def dispositionTask(clip):
    """Decide whether an event is a planet candidate or not

    TODO:
    Much of this should be parcelled off into a function
    """

    #Thresholds
    snrThreshold = clip['config.minSnrForDetection']
    lppThreshold = clip['config.maxLppForTransit']
    minProbForFail = clip['config.minProbDiffImgCentroidForFail']
    snr = clip['trapFit.snr']
    modshiftDict = clip['modshift']
    centroidArray = clip['diffImg.centroid_timeseries']

    out = clipboard.Clipboard(isSignificantEvent=True, isCandidate=True, \
        reasonForFail="None")


    #Compute centroid offset and significance
    centVet = {'Warning':"None"}
    try:
        prob, chisq = cent.measureOffsetProbabilityInTimeseries(centroidArray)
    except ValueError, e:
        centVet['Warning'] = "Probability not computed: %s" %(e)
        prob = 0
        chisq = 0

    centVet['probabilityOfOffset'] = prob
    centVet['chiSquaredOfOffset'] = chisq
    centVet['numCadencesWithCentroids'] = int( np.sum(centroidArray[:,1] > 0))


    centVet['isCentroidFail'] = False
    if np.isfinite(prob):
        if prob > minProbForFail:
            centVet['isCentroidFail'] = True
            out['isCandidate'] = False
            out['reasonForFail'] = "Centroid offset probability is %.1e" %(prob)
    else:
        centVet['Warning'] = "Probability is Nan"

    out['centroidVet'] = centVet


    ####################################################

    #Parse modshift results
    fluxVetDict = RoboVet.roboVet(modshiftDict)
    out['fluxVet'] = fluxVetDict
    assert(fluxVetDict['disp'] in ["candidate", "false positive"])

    if fluxVetDict['disp'] == "false positive":
        out['isCandidate'] = False
        out['reasonForFail'] = fluxVetDict['comments']

        if fluxVetDict['not_trans_like'] > 0:
            out['isSignificantEvent'] = False



    #Check LPP, if it's available
    Tlpp_linear = clip.get('lpp.TLpp', 0)
    if Tlpp_linear > lppThreshold:
        out['isSignificantEvent'] = False
        out['isCandidate'] = False
        out['reasonForFail'] = "TLpp (%.1f) above threshold %.1f" \
            %(Tlpp_linear, lppThreshold)


    #Check SNR
    if snr < snrThreshold:
        out['isSignificantEvent'] = False
        out['isCandidate'] = False
        out['reasonForFail'] = "SNR (%.1f) below threshold %.1f" \
            %(snr, snrThreshold)

    clip['disposition'] = out

    #Enforce contract
    clip['disposition.isSignificantEvent']
    clip['disposition.isCandidate']
    clip['disposition.reasonForFail']
    return clip




import dave.plot.daveplot as daveplot
@task.task
def plotTask(clip):

    time = clip['serve.time']
    raw  = clip['extract.rawLightcurve']
    flux = clip['detrend.flux_frac']
    fl = clip['detrend.flags']
    model = clip['trapFit.bestFitModel']

    epic = clip['value']
    basePath = clip['config.onepageBasename']
    period_days = clip['trapFit.period_days']
    epoch_bkjd = clip['trapFit.epoch_bkjd']
    dur_hrs =  clip['trapFit.duration_hrs']

    basename = getOutputBasename(basePath, epic)
    errorCode = daveplot.onepage(basename, epic, time[~fl], raw[~fl],\
        flux[~fl], model[~fl], period_days, epoch_bkjd, dur_hrs/24.0)

    clip['plot'] = errorCode

    return clip



def saveOnError(clip):
    """Note this is not a task, because it should run even if
    an exception is raised"""

    if 'exception' in clip.keys():
        print "Error found, saving clip..."
        print clip['exception']
        saveClip(clip)
    return clip


import shelve
def saveClip(clip):

    try:
        value = clip['value']
        campaign = clip['config.campaign']
        path = clip.get('config.clipSavePath', ".")
        path = os.path.join(path, "%09i" %value)

        if not os.path.exists(path):
            os.mkdir(path)

        #The problem with this is how do I which tasks to run
        #when I restore?
        keysToSkip = clip.get('config.keysToIgnoreWhenSaving', [])

        fn = "c%09i-%02i-%02i-%04i-%s.clip" %(value,campaign,np.round(clip.bls.period*10),np.round(clip.bls.epoch),clip.config.detrendType)
        fn = os.path.join(path, fn)

        #Remove old clip files and any clip.db, clip.bak etc.
        import glob
        fList = glob.glob("%s*" %(fn))
        for f in fList:
            os.remove(f)

        sh = shelve.open(fn)
        for k in clip.keys():
            if k in keysToSkip:
                sh[k] = "Clip not saved"
            else:
                sh[k] = clip[k]
        sh.close()
    
    except Exception, e:
        print "WARN: Error in saveClip: %s" %(e)
        clip['exception'] = str(e)

    return clip

def loadTpfAndLc(k2id, campaign, ar, detrendType):
#    ar = mastio.K2Archive(storeDir)  #Removed by SEM to generalize this function

    out = dict()
    fits, hdr = ar.getLongTpf(k2id, campaign, header=True, mmap=False)
    hdr0 = ar.getLongTpf(k2id, campaign, ext=0, mmap=False)
    cube = tpf.getTargetPixelArrayFromFits(fits, hdr)

    out['cube'] = cube
    out['tpfHeader'] = hdr
    out['tpfHeader0'] = hdr0

    fits, hdr2 = ar.getLongCadence(k2id, campaign, header=True)
    data = kplrfits.getNumpyArrayFromFitsRec(fits)
    lookup = """ TIME TIMECORR CADENCENO
                 SAP_FLUX SAP_FLUX_ERR SAP_BKG SAP_BKG_ERR
                 PDCSAP_FLUX PDCSAP_FLUX_ERR SAP_QUALITY
                 PSF_CENTR1 PSF_CENTR1_ERR PSF_CENTR2 PSF_CENTR2_ERR
                 MOM_CENTR1 MOM_CENTR1_ERR MOM_CENTR2 MOM_CENTR2_ERR
                 POS_CORR1 POS_CORR2""".split()
    data = nca.Nca(data)
    data.setLookup(1, lookup)


    #Load lightcurves from a specific detrending, and replace
    #the pdc time series with the new detrending
    key = detrendType.upper()
    if key == "PDC":
        pass

    elif key == "EVEREST":
        detrendAr = mastio.EverestArchive()
        fits2 = detrendAr.getLongCadence(k2id, campaign)
        flux = fits2['FLUX']
        assert len(flux) == len(data)
        data[:, 'PDCSAP_FLUX'] = flux

    elif key == "SFF":
        detrendAr = mastio.VanderburgArchive()
        lightcurve = detrendAr.getLongCadence(k2id, campaign)
        sffTime = lightcurve[:,0]
        flux = lightcurve[:, 1]

        idx = mapTime2ToIndexOfTime1(data[:, 'TIME'], sffTime)
        data[idx, 'PDCSAP_FLUX'] = flux
        data[~idx, 'PDCSAP_FLUX'] = np.nan

    elif key == "AGP":
        detrendAr = mastio.K2SCArchive()
        agpFits = detrendAr.getLongCadence(k2id, campaign)
        agpTime = agpFits['time']
        agpFlux = agpFits['flux']

        idx = mapTime2ToIndexOfTime1(data[:, 'TIME'], agpTime)
        data[idx, 'PDCSAP_FLUX'] = agpFlux
        data[~idx, 'PDCSAP_FLUX'] = np.nan
    else:
        raise ValueError("Unrecognised detrending type %s" %(detrendType))


    out['time'] = fits['TIME'].copy()
    out['flags'] = fits['SAP_QUALITY'].copy()
    out['socData'] = data
    return out


def mapTime2ToIndexOfTime1(time1, time2):
    """Compute the indices of time1 that correspond closest to values of time2

    Some K2 detrended lightcurves (like PDC) contain elements for every
    cadence, others (like sff) only contain elements for good data.
    This function enables you to figure out which rows of time1 correspond
    to the same times in time2, so you can place place the first and second
    data sets in contiguous arrays

    Typical Usage:
    ----------------
    A typical usage would look like
    ::
        idx = mapTime2ToIndexOfTime1(data1[:, 'TIME'], data2[:, 'TIME'])
        data1[idx, 'FLUX2'] = data2[:, 'FLUX']
        data1[~idx, 'FLUX2'] = BAD_VALUE

    Inputs:
    ------------
    time1:
        (1d np array) Array of times to map into
    time2
        (1d np array) Array of times to map from. Typically ``len(time2) < len(time1)``

    Returns:
    -----------
    A boolean array of length ``time1``

   """

    t1 = np.atleast_2d(time1)
    t2 = np.atleast_2d(time2)
    dt = t1 - t2.transpose()
    dt = np.fabs(dt.asarray())
    assert dt.ndim == 2
    dt = np.nanmin(dt, axis=0)
    idx = dt < 1e-8  #Two times must agree very well to be accepted
    assert len(idx) == len(time1)
    return idx



def getOutputBasename(basePath, epic):
    """Get the output basename for any files a task creates

    Inputs:
    ----------
    basePath
        (string) Path where all files should be output
    epic
        (int or string) Epic number of star


    Returns:
    -----------
    (string) a basename for files.


    Outputs:
    ----------
    Function attempts to create output directory if it doesn't exist.

    Example:
    -------------
    >>> getOutputBasename("/home/dave/c6", 123456789)
    "/home/dave/c6/123456789/123456789"

    The calling task can then create a bunch files like
    "/home/dave/c6/123456789/123456789-file1.txt"
    "/home/dave/c6/123456789/123456789-image1.png" etc.
    """

    epicStr = str(int(epic))
    path = os.path.join(basePath, epicStr)

    if not os.path.exists(path):
        os.mkdir(path)
    if not (os.access(path, os.W_OK)):
        raise IOError("Can't write to output directory %s" %path)

    return os.path.join(path, epicStr)


@task.task
def cotrendSffDataTask(clip):
    """Produce a cotrended lightcurve in units of fractional amplitude"""
    from dave.detrendThis.detrendThis import detrendThat

    time = clip['serve.time']
    flags = clip['extract.flags']
    xbar = clip['extract.centroid_col']
    ybar = clip['extract.centroid_row']
    rawflux = clip['extract.rawLightcurve']

    outcorflux, outcorflatflux, outcorrection, detrendFlags = detrendThat(
                time[~flags], rawflux[~flags], xbar[~flags], ybar[~flags],
                ferr=None,
                qflags=None,
                inpflag=None,
                ap=4.0)

    newDetrendFlags = np.zeros_like(flags)
    newDetrendFlags[~flags] = ~detrendFlags
    flags |= newDetrendFlags
    flux = rawflux.copy()
    flux[~flags] = outcorflux

    newCorr = np.zeros_like(rawflux)
    newCorr[~flags] = outcorrection

    #Cotrending may also produce Nans
    flags |= ~np.isfinite(flux)

    #Remove dc offset
    dcOffset = np.median( flux[~flags])
    flux = (flux/ dcOffset) - 1
    clip['cotrend'] = {'flux_frac': flux}
    clip['cotrend.dcOffset'] = dcOffset
    clip['cotrend.flags'] = flags
    clip['cotrend.dcOffset'] = dcOffset
    clip['cotrend.source'] = "SFF Cotrend"
    clip['cotrend.correction'] = newCorr

    #Enforce contract
    clip['cotrend.flux_frac']
    return clip


@task.task
def extractLightcurveFromTpfTask(clip):
    from dave.tpf2lc.tpf2lc import optimalAperture
    time = clip['serve.time']
    fluxcube = clip['serve.cube']
    data = clip['serve.socData']
    socflux = data[:,'SAP_FLUX']
    numInitialCadencesToIgnore = clip['config.numInitialCadencesToIgnore']
    flagValues = clip.get('serve.flags', data[:, 'SAP_QUALITY'])

    #Convert flags to a boolean, and flag other bad data
    mask = kplrfits.getMaskForBadK2Data()
    flags = (flagValues & mask).astype(bool)
    flags |= ~np.isfinite(time)
    flags |= ~np.isfinite(socflux)
    flags[socflux<1] = True
    flags[:numInitialCadencesToIgnore] = True

    newtime, flux, xbar, ybar, _, _ = optimalAperture(time[~flags], fluxcube[~flags], flags[~flags],
        qual_cut=False,
        bg_cut=4)

    newY = socflux.copy()
    newXbar = np.zeros_like(socflux)
    newYbar = np.zeros_like(socflux)
    newY[~flags] = flux
    newXbar[~flags] = xbar
    newYbar[~flags] = ybar

    #Placeholder. Use the SOC PA data for the lightcurve
    out = dict()
    out['rawLightcurve'] = newY
    clip['extract'] = out
    clip['extract.source'] = "Labeled Extraction Pipeline"
    clip['extract.flags'] = flags
    clip['extract.centroid_col'] = newXbar
    clip['extract.centroid_row'] = newYbar

    #Enforce contract
    clip['extract.rawLightcurve']
    clip['extract.flags']
    return clip

import dave.milesplay.productionPCA2 as milesPCA
@task.task
def getMilesLightcurveTask(clip):
    
    # get input data from clip
    cube = clip['serve.cube']
    time = clip['serve.time']
    flags = clip['serve.flags'] > 0
    flags |= ~np.isfinite(time)

    nt, nr, nc = cube.shape
    pixSeries, pixSeriesNorm = milesPCA.pixSeriesPrep(cube, flags)
    
    totLightcurve = milesPCA.getRawLightcurve(pixSeries)
    
    # package results
    out = clipboard.Clipboard()
    out['rawLightcurve'] = totLightcurve
    out['time'] = time
    out['pixSeries'] = pixSeries
    out['pixSeriesNorm'] = pixSeriesNorm
    out['flags'] = flags
    
    clip['extract'] = out
    
    #Enforce contract
    clip['extract.rawLightcurve']
    clip['extract.time']
    clip['extract.pixSeries']
    clip['extract.pixSeriesNorm']
    clip['extract.flags']
    return clip
    
@task.task    
def milesCotrendDataTask(clip):
    
    #pixSeries = clip['extract.pixSeries']
    pixSeriesNorm = clip['extract.pixSeriesNorm']
    rawLightcurve = clip['extract.rawLightcurve']
    flags = clip['extract.flags']
    numPixels = pixSeriesNorm.shape[0]
    evalues, prinComps = milesPCA.performSVD(pixSeriesNorm)
    optimalNumPC = milesPCA.chooseNumPrinComps(prinComps, rawLightcurve, numPixels, evalues)
    bestFit = milesPCA.curveFit(prinComps[:optimalNumPC], rawLightcurve)
    print  optimalNumPC
    
    bestCorrectedLightcurve = rawLightcurve - bestFit
    flux_frac = bestCorrectedLightcurve
    #flux_frac = flux_frac - np.mean(flux_frac)
    flux = flux_frac/np.mean(flux_frac) - 1
    flux -= np.mean(flux)
    #flux = flux[~flags]
    #Remove dc offset
    #dcOffset = np.mean( flux_frac)
    #flux = (flux_frac/ dcOffset) - 1
    clip['cotrend'] = {'flux_frac': flux}
    #clip['cotrend.dcOffset'] = dcOffset
    clip['cotrend.flags'] = flags
    clip['cotrend.source'] = "Pixel Level PCA"
    clip['cotrend.correctedCurve'] = bestCorrectedLightcurve
    clip['cotrend.numPC'] = optimalNumPC
    
    
    # Enforce Contract
    clip['cotrend.flux_frac']
    clip['cotrend.correctedCurve']
    clip['cotrend.numPC']
    
    return clip