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
    cfg['nPointsForMedianSmooth'] = 2*48
    cfg['blsMinPeriod'] = 0.5
    cfg['blsMaxPeriod'] = 30

    #Vetting parameters
    cfg['minSnrForDetection'] = 10.
    cfg['maxLppForTransit'] = 0.03
    #How significant must centroid offset be to claim a false positive.
    #This value is between 0 and 1. Higher values mean fewer false positives
    cfg['minProbDiffImgCentroidForFail'] = 0.99


    #The LPP mapping file is used by LPP to define the regions of param
    #space where the transits cluster.
    path = lpp.getLppDir()
    cfg['lppMapFilePath'] = os.path.join(path, "octave/maps/mapQ1Q17DR24-DVMed6084.mat")

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
        computeCentroidsTask rollPhaseTask cotrendDataTask detrendDataTask
        blsTask trapezoidFitTask vetTask plotTask""".split()   # Jeff added plotTask
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

    clip['serve'] = loadTpfAndLc(k2id, campaign, storeDir)

    #Enforce contract. (Make sure expected keys are in place)
    clip['serve.time']
    clip['serve.cube']
    clip['serve.socData']
    clip['serve.tpfHeader']
    return clip


@task.task
def extractLightcurveTask(clip):
    data = clip['serve.socData']

    flagValues = clip.get('serve.flags', data[:, 'SAP_QUALITY'])
    flux = data[:, 'SAP_FLUX']

    #Convert flags to a boolean.
    mask = kplrfits.getMaskForBadK2Data()
    flags = (flagValues & mask).astype(bool)

    #Flag bad values
    flags[~np.isfinite(flux)] = True
    flags[flux<1] = True

    #Placeholder. Use the SOC PA data for the lightcurve
    out = dict()
    out['rawLightcurve'] = flux
    clip['extract'] = out
    clip['extract.source'] = "SOC PA Pipeline"
    clip['extract.flags'] = flags

    #Enforce contract
    clip['extract.rawLightcurve']
    return clip


@task.task
def cotrendDataTask(clip):
    """Produce a cotrended lightcurve in units of fractional amplitude"""

    data = clip['serve.socData']
    flags = clip['extract.flags']
    time = data[:, 'TIME']
    flux = data[:, 'PDCSAP_FLUX']

    flags |= ~np.isfinite(time)
    flags |= ~np.isfinite(flux)

    #Remove dc offset
    dcOffset = np.median( flux[~flags])
    flux = (flux/ dcOffset) - 1
    clip['cotrend'] = {'flux_frac': flux}
    clip['cotrend.dcOffset'] = dcOffset
    clip['cotrend.flags'] = flags
    clip['cotrend.source'] = "SOC PDC Pipeline"

    #Enforce contract
    clip['cotrend.flux_frac']
    return clip


@task.task
def detrendDataTask(clip):
    flux = clip['cotrend.flux_frac']
    flags = clip['cotrend.flags']

    nPoints = clip['config.nPointsForMedianSmooth']

    #When you detrend, you must do something about the gaps and bad values.
    #This is the simplest possible thing. Replace all bad/missing data with
    #zeros. This is a placehold. Bad data inside a transit is replaced with
    #a zero, which is not what you want.
    flux[flags] = 0

    #Do a simple detrend.
    detrend = kplrfits.medianSubtract1d(flux, nPoints)
    clip['detrend'] = dict()
    clip['detrend.flux_frac'] = detrend
    clip['detrend.flags'] = flags
    clip['detrend.source'] = "Simple Median detrend"

    assert(detrend is not None)
    return clip


@task.task
def computeCentroidsTask(clip):
    data = clip['serve.socData']

    cent_colrow = np.empty( (len(data), 2))
    cent_colrow[:,0] = data[:, 'MOM_CENTR1']
    cent_colrow[:,1] = data[:, 'MOM_CENTR2']
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


import dave.blsCode.bls_ktwo as bls
@task.task
def blsTask(clip):
    time_days = clip['serve.time']
    flux_norm = clip['detrend.flux_frac']
    flags = clip['detrend.flags']
    minPeriod = clip['config.blsMinPeriod']
    maxPeriod = clip['config.blsMaxPeriod']

    #Zero out the bad data. This crashes BLS
#    flux_norm[flags] = 0
#    assert(np.all( np.isfinite(flux_norm)))

    idx = flags == 0
    period, epoch, duration, depth, bls_search_periods, convolved_bls = \
        bls.doSearch(time_days[idx], 1+flux_norm[idx], minPeriod, maxPeriod)

    out = clipboard.Clipboard()
    out['period'] = period
    out['epoch'] = epoch
    out['duration_hrs'] = duration * 24
    out['depth'] = depth
    out['bls_search_periods'] = bls_search_periods
    out['convolved_bls'] = convolved_bls
    clip['bls'] = out

    #Enforce contract
    clip['bls.period']
    clip['bls.epoch']
    clip['bls.duration_hrs']
    return clip


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

    #Place holder, use Susan's version when it shows up.
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
    depth_frac = clip['bls.depth']

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


import dave.trapezoidFit.trapfit as trapFit
import dave.vetting.ModShift as ModShift
@task.task
def modshiftTask(clip):

    time = clip['serve.time']
    flux = clip['detrend.flux_frac']
    fl = clip['detrend.flags']

    epic = clip['value']
    basename = clip['config.modshiftBasename'] #  Jeff edited to remove ->  + "%010i" %(epic)   since modshift cpp code already appends basename
    basename = "%s%010i" %(basename, epic)  # Jeff modified

    period_days = clip['trapFit.period_days']
    epoch_bkjd = clip['trapFit.epoch_bkjd']
    dur_hrs =  clip['trapFit.duration_hrs']
    ingress_hrs = clip['trapFit.ingress_hrs']
    depth_ppm = 1e6*clip['trapFit.depth_frac']
    objectname = "EPIC " + str(epic)  # Name that will go in title of modshift plot

    subSampleN= 15
    ioBlock = trapFit.trapezoid_model_onemodel(time[~fl], period_days, \
        epoch_bkjd, depth_ppm, dur_hrs, \
        ingress_hrs, subSampleN)
    model = ioBlock.modellc -1   #Want mean of zero
    #model *= -1  #Invert for testing

    out = ModShift.runModShift(time[~fl], flux[~fl], model, basename, \
        objectname, period_days, epoch_bkjd)

    clip['modshift'] = out

    #I don't know which values are important, so I can't enfornce contract yet
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



    #Check LPP, if it's available
    Tlpp_linear = clip.get('lpp.TLpp', 0)
    if Tlpp_linear > lppThreshold:
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

    epic = clip['value']
    basename = clip['config.onepageBasename'] #  Jeff removed end so it isn't redudant with plotting module epic adding + "%010i" %(epic)
    period_days = clip['trapFit.period_days']
    epoch_bkjd = clip['trapFit.epoch_bkjd']
    dur_hrs =  clip['trapFit.duration_hrs']
    ingress_hrs = clip['trapFit.ingress_hrs']
    depth_ppm = 1e6*clip['trapFit.depth_frac']

    basename = "%s%010i" %(basename, epic)

    subSampleN= 15
    ioBlock = trapFit.trapezoid_model_onemodel(time[~fl], period_days, \
        epoch_bkjd, depth_ppm, dur_hrs, \
        ingress_hrs, subSampleN)
    model = ioBlock.modellc -1   #Want mean of zero

    out = daveplot.onepage(basename,epic,time[~fl],raw[~fl],flux[~fl],model,period_days,epoch_bkjd,dur_hrs/24.0)

    clip['plot'] = out

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

        #The problem with this is how do I which tasks to run
        #when I restore?
        keysToSkip = clip.get('config.keysToIgnoreWhenSaving', [])

        fn = "c%09i-%02i.clip" %(value, campaign)
        fn = os.path.join(path, fn)

        if os.path.exists(fn):
            os.remove(fn)

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

def loadTpfAndLc(k2id, campaign, storeDir):
    ar = mastio.K2Archive(storeDir)

    out = dict()
    fits, hdr = ar.getLongTpf(k2id, campaign, header=True)
    hdr0 = ar.getLongTpf(k2id, campaign, ext=0)
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
    out['socData'] = data
    out['time'] = fits['TIME']
    out['flags'] = fits['SAP_QUALITY']
    return out




