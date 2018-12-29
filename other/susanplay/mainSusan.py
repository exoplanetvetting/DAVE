
"""
This is a template top level script.

Please don't edit this file. Instead, copy it to
youname_main.py, then run and edit that file.

"""

import dave.pipeline.pipeline as dpp
import dave.pipeline.clipboard as clipboard

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
    cfg['debug'] = False


    tasks = """dpp.checkDirExistTask dpp.serveTask dpp.extractLightcurveFromTpfTask
        dpp.computeCentroidsTask dpp.rollPhaseTask dpp.cotrendSffDataTask
        dpp.detrendDataTask dpp.fblsTask dpp.trapezoidFitTask dpp.lppMetricTask dpp.modshiftTask
        dpp.measureDiffImgCentroidsTask dpp.dispositionTask
        dpp.saveClip""".split()   
    cfg['taskList'] = tasks
    
    searchTaskList = """blsTask trapezoidFitTask modshiftTask
    measureDiffImgCentroidsTask dispositionTask""".split()

    cfg['searchTaskList'] = searchTaskList

    return cfg



def runOne(k2id, config):
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
    for t in taskList:
        f = eval(t)

    #Now run them.
    for t in taskList:
        f = eval(t)
        clip = f(clip)

    return clip



def runOneEphem(k2id,period,epoch,config,duration=2,depth=.0001):
    """
    Run just the vetting and return an output.
    Inputs:
    -------------
    k2id
        (int) Epic id of the target to run on.
    period
        (float) period of the target
    epoch
        (float) Time in days
    config
        (dict) Dictionary of configuration parameters
       
    """
    
    tasks = """dpp.checkDirExistTask dpp.serveTask dpp.extractLightcurveFromTpfTask
        dpp.computeCentroidsTask dpp.rollPhaseTask dpp.cotrendSffDataTask
        dpp.detrendDataTask dpp.trapezoidFitTask dpp.lppMetricTask 
        dpp.modshiftTask dpp.measureDiffImgCentroidsTask dpp.dispositionTask
        dpp.saveOnError""".split()  
    
    taskList = tasks;

    clip = clipboard.Clipboard()
    clip['config'] = config
    clip['value'] = k2id
    out = clipboard.Clipboard()
    out['period'] = period
    out['epoch'] = epoch
    out['duration_hrs'] = duration
    out['depth'] = depth
    clip['bls'] = out
            
    

    #Check that all the tasks are properly defined
    for t in taskList:
        f = eval(t)

    #Now run them.
    for t in taskList:
        f = eval(t)
        clip = f(clip)

    return clip