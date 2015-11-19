# -*- coding: utf-8 -*-


from multiprocessing import pool
import multiprocessing
import contextlib

import clipboard
import traceback
import parmap
import time
import pdb
import sys

__version__ = "$Id: task.py 2130 2015-09-11 16:56:55Z fmullall $"
__URL__ = "$URL: svn+ssh://flux/home/fmullall/svn/kepler/k2phot/task.py $"



"""
A task is a serial element of a pipeline intended to be run in parallal.
It contains a specific step of the process and helps make your pipeline
more modular. It also makes debugging error easier.

Here is an example

@task
def exampleTask(clip):
    #clip is a Clipboard object
    cache = clip['config.cachePath']
    outPath = clip['config.outputPath']

    clip['example'] = runSillyOperation(cache, out)

    #Check required key is produced
    clip['example.dummyValue']
    return clip

def runSillyOperation(cache, out):
    out = dict()
    out['dummyValue'] = 42
    return out
"""



def task(func):
    """A decorator for a task function

    The k2phot model is that you write wrapper code around your pipeline functions
    that extract values from a clipboard to pass as input arguments,
    and store the return values back in the clipboard. This wrapper code
    is called a "task".

    This decorator watches for any exceptions thrown by the task (or
    underlying code) and decides whether to fire up the debugger to
    figure out what's going on, or to fail gracefully by merely storing
    the raised exception for later debugging. This decorator considerably reduces
    the code duplication between different tasks.

    To use this code in your function, use the following import statement
    from task import task
    then simply write the text "@text" above your task. For example

    from task import task

    @task
    def exampleTask(clip):
        pass

    See an example down below.
    """

    def wrapperForTask(*args, **kwargs):
        """Decorator for a k2phot task

        Catches exceptions and either stores them or passes them to the debugger
        """
        assert(len(args) == 1)
        assert(len(kwargs) == 0)

        clip = args[0]
        if 'exception' in clip.keys():
            print "INFO: %s not run because exception previously raised" \
                %(func.func_name)
            return clip

        if "__meta__" not in clip.keys():
            clip["__meta__"] = dict()

        debug = clip.get('config.debug', False)
        t0 = time.time()
        try:
            clip = func(*args, **kwargs)
        except SyntaxError, e:
            raise(e)
        except Exception, e:
            if debug:
                print e
                pdb.post_mortem(sys.exc_info()[2])
                raise e
            else:
                clip['exception'] = e
                clip['backtrace'] = traceback.format_exc()

        if not isinstance(clip, clipboard.Clipboard):
            if not isinstance(clip, dict):
                throwable = ValueError("FAIL: %s did not return a clipboard" %(func.func_name))

                if debug:
                    raise throwable
                else:
                    print throwable
                    clip = {'exception': throwable}
                    clip['functionName'] = func.__name__

        key = "%s-elapsedTime" %(func.__name__)
        clip["__meta__"][key] = time.time() - t0
        return clip

    wrapperForTask.__name__ = func.__name__
    return wrapperForTask



def runAll(func, iterable, config):
    """Run func over every element on iterable in parallel.

    Not yet run or tested.

    Inputs:
    ----------
    func
	(A function) The top level function, e.g runOne(), below

    iterable
	(list, array, etc.) A list of values to operate on.

    config
	(Clipboard) A configuration clipboard.
    """

    count = multiprocessing.cpu_count() - 1
    p = pool.Pool(count)


    parallel = config.get('debug', False)

    with contextlib.closing(pool.Pool(count)) as p:
        out = parmap.map(runOne, iterable, config, pool=p, parallel=parallel)

    return out




def runOne(value, config):
    """A sample top level function to run a single instance of a pipeline

    Not yet run or tested.

    Inputs:
    ----------
    value
	(undefined)  The input value that is unique to this run, eg the name of
	the file to process, or the value to perform some computation on.

    config
	(Clipboard) A configuration clipboard. Must contain the key 'taskList'. See below


    Returns:
    -----------
    A clipboard, containing the results of the processing.

    Notes:
    ----------
    The config clipboard must contain the key taskList. taskList is
    a list of strings, each listing the function name of a task to
    run in the pipeline, in the order they should be run. The
    first task should look of the unique value in clip['value']

    The tasks as passed as strings because functions can't bepa
    ssed to parallel processes.This means that your tasks must
    be defined in scope of this file. The easiest way
    to do this looks like

    import task
    import stuff

    task.func1 = stuff.func1
    task.func2 = stuff.func2

    config['taskList'] = ["func1", "func2"]
    """


    taskList = config['taskList']

    clip = clipboard.Clipboard()
    clip['config'] = config
    clip['value'] = value

    for t in taskList:
        f = eval(t)
        clip = f(clip)

    return clip
