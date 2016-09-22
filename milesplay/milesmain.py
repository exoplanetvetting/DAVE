
"""
This is a template top level script.


"""

import dave.pipeline.pipeline as dpp
import dave.pipeline.clipboard as clipboard
import dave.pipeline.task as task
import numpy as np
import gc

def main():
    """A bare bones main program"""
    cfg = loadMyConfiguration()
#==============================================================================
#     campaignList = [3,5,3,3,3,3,3,3,3,3]
#     epicList = [206103150, 211351816, 206348688, 206268299, 206247743, 206245553
#                 ,206192813,206181769 , 206162305 , 206159027    ]
#     
#     candInfo = []
#     with open('/Users/Miles/seti/candidates_to_run.tsv') as file:
#         for line in file:
#             line = line.strip("\n")
#             line = line.split("\t")
#             line[1]  = line[1].strip("EPIC ")
#             candInfo.append(line)
#     candInfo = np.array(candInfo)
#     epicList = candInfo.T[1][1:].astype(int)
#     campaignList = candInfo.T[2][1:].astype(int)
#==============================================================================
    file = '/Users/Miles/seti/c8_stars.csv'
    with open(file) as f:
        infoMatrix = []
        for l in f:
            l = l.strip('\n')
            l = l.split(',')
            infoMatrix.append(l)
    infoMatrix = np.array(infoMatrix)[1:]
    epics_mags = np.column_stack((infoMatrix.T[0,:].astype(int), infoMatrix.T[3,:].astype(float)))    
    sorted_epicsMags =  epics_mags[epics_mags[:,1].argsort()]
    epicList = sorted_epicsMags.T[0]
    campaignList = np.ones_like(epicList)*8
    n=0
    for epic, campaign in zip(epicList[:10], campaignList[:10]):
        print '\nCOUNTER:',n
        print "EPIC: %i CAMPAIGN: %i"%(epic, campaign)
        cfg['campaign'] = campaign
        runOne(epic, cfg)
        n+=1


def loadMyConfiguration():
    """Load the default pipeline configuration and adjust as necessary
    """

    cfg = dpp.loadDefaultConfig()

    cfg['lppMapFilePath'] = '/Users/Miles/seti/dave/lpp/octave/maps/mapQ1Q17DR24-DVMed6084.mat'

    cfg['detrendType'] = 'PDC'    
    #Edit the default configuration to your taste.
    #Change anything else you don't like about the default config here.
    cfg['debug'] = False
    cfg['timeout_sec'] = 0

    tasks = """dpp.checkDirExistTask dpp.serveTask dpp.getMilesLightcurveTask
        dpp.computeCentroidsTask dpp.rollPhaseTask dpp.milesCotrendDataTask
        dpp.detrendDataTask dpp.fblsTask dpp.trapezoidFitTask  dpp.saveClip""".split()
        
#==============================================================================
#         dpp.modshiftTask
#         dpp.lppMetricTask dpp.measureDiffImgCentroidsTask dpp.dispositionTask
#         
#==============================================================================
        
    #tasks = ["dpp.getMilesLightcurveTask"]
    cfg['taskList'] = tasks
    cfg['numInitialCadencesToIgnore'] = 100

#    searchTaskList = """placeholderBls trapezoidFitTask modshiftTask
#    measureDiffImgCentroidsTask dispositionTask""".split()
    searchTaskList = """blsTask trapezoidFitTask modshiftTask
    measureDiffImgCentroidsTask dispositionTask""".split()

    cfg['searchTaskList'] = searchTaskList
    return cfg


import multiprocessing
import contextlib
import dave.pipeline.parmap as parmap
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

main()
