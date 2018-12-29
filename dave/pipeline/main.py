
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
    cfg['debug'] = True

    tasks = """dpp.serveTask dpp.extractLightcurveTask
        dpp.computeCentroidsTask dpp.rollPhaseTask dpp.cotrendDataTask
        dpp.detrendDataTask dpp.fblsTask dpp.trapezoidFitTask
        dpp.modshiftTask dpp.measureDiffImgCentroidsTask dpp.dispositionTask
        dpp.plotTask dpp.saveOnError""".split()

#    tasks = """dpp.checkDirExistTask dpp.serveTask dpp.extractLightcurveTask
#        dpp.computeCentroidsTask dpp.rollPhaseTask dpp.cotrendDataTask
#        dpp.detrendDataTask dpp.blsTask dpp.trapezoidFitTask dpp.dispositionTask
#        dpp.plotTask dpp.saveOnError""".split()   # Jeff added
    cfg['taskList'] = tasks

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


