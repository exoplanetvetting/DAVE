
"""
This is a template top level script.


"""

import dave.pipeline.pipeline as dpp
import dave.pipeline.clipboard as clipboard

import gc

def main():
    """A bare bones main program"""
    cfg = loadMyConfiguration()

    epicList = [206103150]

    for epic in epicList:
        runOne(epic, cfg)


def loadMyConfiguration():
    """Load the default pipeline configuration and adjust as necessary
    """

    cfg = dpp.loadDefaultConfig()

    #Edit the default configuration to your taste.
    #Change anything else you don't like about the default config here.
    cfg['debug'] = True
    cfg['timeout_sec'] = 0

    tasks = """dpp.checkDirExistTask dpp.serveTask dpp.extractLightcurveTask
        dpp.computeCentroidsTask dpp.rollPhaseTask dpp.cotrendDataTask
        dpp.detrendDataTask dpp.fblsTask dpp.trapezoidFitTask dpp.modshiftTask
        dpp.lppMetricTask dpp.measureDiffImgCentroidsTask dpp.dispositionTask
        dpp.saveClip""".split()


    cfg['taskList'] = tasks


#    searchTaskList = """placeholderBls trapezoidFitTask modshiftTask
#    measureDiffImgCentroidsTask dispositionTask""".split()
    searchTaskList = """blsTask trapezoidFitTask modshiftTask
    measureDiffImgCentroidsTask dispositionTask""".split()

    cfg['searchTaskList'] = searchTaskList
    return cfg


import multiprocessing
import contextlib
import parmap
from multiprocessing import pool
def runAll(func, iterable, config):
    """Run func over every element on iterable in parallel.

    Inputs:
    ----------
    func
	(A function) The top level function, e.g runOne(), below

    iterable
	(list, array, etc.) A list of values to operate on.

    config
	(Clipboard) A configuration clipboard.
    """

    #Turn off parallel when debugging.
    parallel = True
    if config.get('debug', False):
        parallel = False

    count = multiprocessing.cpu_count()-1

    with contextlib.closing(pool.Pool(count)) as p:
        out = parmap.map(func, iterable, config, pool=p, parallel=parallel, chunksize=5)

    return out


def runOne(k2id, config, returnClip=False):
    """Run the pipeline on a single target.

    Inputs:
    ------------
    k2id
        (int) Epic id of target to run on.

    config
        (dict) Dictionary of configuration parameters as created by, e.g
        loadMyConfiguration()

    Returns:
    ---------
    A clipboard containing the results.

    Notes:
    ---------
    Don't edit this function. The pipeline can recover gracefully from
    errors in any individual task, but an error in this function will crash
    the pipeline
    """

    taskList = config['taskList']

    clip = clipboard.Clipboard()
    clip['config'] = config
    clip['value'] = k2id

    #Check that all the tasks are properly defined
    print "Checking tasks exist"
    for t in taskList:
        f = eval(t)

    #Now run them.
    for t in taskList:
        print "Running %s" %(t)
        f = eval(t)
        clip = f(clip)

    gc.collect()
    if returnClip:
        return clip




import dave.misc.noise as noise
import dave.fileio.kplrfits as kplrfits
import task
@task.task
def newDetrendDataTask(clip):
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

    rollTweakAmp = noise.computeRollTweakAmplitude(detrend[~flags])
    clip['detrend.rollTweakAmp'] = rollTweakAmp

    cdpp6 = noise.computeSgCdpp_ppm(detrend[~flags])
    clip['detrend.cdpp6_ppm'] = cdpp6

    perPointScatter = noise.estimateScatterWithMarshallMethod(detrend[~flags])
    clip['detrend.perPointScatter_ppm'] = 1e6*perPointScatter


    assert(detrend is not None)
    return clip
