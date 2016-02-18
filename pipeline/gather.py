# -*- coding: utf-8 -*-
"""
Utility functions to sweep through a list of saved clipboards
and *gather* information, or compute quantities. The two functions
of interest are:

gatherValues:  Extract a value stored in each clipboard
gatherFunction Call a function on each clipboard.


Examples:
clipList = glob.glob("clips/*.clip")
(epics, periods) = gatherValues(clipList, 'trapFit.period_days')


def planetStarRadiusRatio(clip):
    return np.sqrt(clip['trapFit.depth_frac'])

(epics, RpRs) = sweepFunction(clipList, planetStarRadiusRatio)
"""


import dave.pipeline.clipboard as dpc
#from multiprocessing import pool
#import multiprocessing
#import contextlib
#import parmap


def gatherValue(clipList, key):
    """Retreive a value from a list of saved clipboard files

    Inputs:
    -----------
    clipList
        A list of filenames of clipboards to read
    key
        (str) The key to search

    Returns:
    ----------
    A tuple of two lists.
    * A list of epics for the clipboards for which values were retreived
    * A list of the values of the input key for each clip

    Notes:
    ---------
    This is a thin wrapper for sweepFunction()

    """

    func = lambda x: x[key]
    return gatherFunction(clipList, func)


#Of no use
#def sweepFunctionA(clipList, func):
#    out = []
#    for c in clipList:
#        clip = dpc.loadClipboard(c)
#        out.append(func(clip))
#    return out


def gatherFunctionB(clipList, func, applyFilter=True):
    """Retreive a value from a list of saved clipboard files

    Inputs:
    -----------
    clipList
        A list of filenames of clipboards to read

    func
        The function to apply to each clipboard. The function should
        have the signature func(clip). It must accept one and only one
        argument.

    Optional Inputs:
    -----------------
    applyFilter
        If True, remove from output any clips for which a value was not
        computed (e.g because the clipboard couldn't be loaded). If False,
        failed results are represented by **None**

    Returns:
    ----------
    A tuple of two lists.

    * A list of epics for the clipboards for which values were retreived
    * A list of the return values of the input function for each clipboard

    Notes:
    -----------
    KeyErrors in func() are caught and handled. All other errors
    raise an exception. I may be more tolerant of other errors in future,
    but right now I want to stop computation if func crashes.
    """
    result = map(lambda x: runFuncOnClip(x, func), clipList)

    #Strip out the values that failed?
    if applyFilter:
        result = filter(lambda x: x[1] is not None, result)

    #Split out the epics and the results
    values = map(lambda x: x[0], result)
    result = map(lambda x: x[1], result)
    return values, result


#def gatherFunctionC(clipList, func, parallel=True, count=None):
#    """This doesn't work. Pickling a function is not allowed.
#
#    I could make it work with some eval hacking. eg
#    import sweep
#    sweep.funcToEval = mySpecialFunc
#    sweep.sweepFunctionC(clipList)
#
#    Then runFuncOnClip(clipFilename) would call
#    funcToEval(clip), instead of being passed in the function.
#    """
#    raise NotImplementedError("This func not implemented. Try sweepFunctionB")
#    if count is None:
#        count = multiprocessing.cpu_count()-1
#
#    with contextlib.closing(pool.Pool(count)) as p:
#        out = parmap.map(runFuncOnClip, clipList, func, pool=p, parallel=parallel)
#
#    return out


gatherFunction = gatherFunctionB

def runFuncOnClip(clipFilename, func):
    """Load a clipboard and apply a function to it

    Inputs:
    ----------
    clipFilename
        (string) Name of clipboard to load
    func
        (function) Function to apply. func must take one argument, a clipboard

    Returns:
    ------------
    (value, result). Value is clip['value'], result = func(clip)
    result is set to None for some errors, other errors in func() raise
    an exception
    """

    try:
        clip = dpc.loadClipboard(clipFilename)
    except IOError, e:
        print e
        return 0, None

    try:
        value = clip['value']
        result = func(clip)
    except KeyError, e:
        print "KeyError on %s: %s" %(clipFilename, e)
        return value, None

    return value, result