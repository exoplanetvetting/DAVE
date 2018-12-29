# Copyright 2017-2018 Orbital Insight Inc., all rights reserved.
# Contains confidential and trade secret information.
# Government Users: Commercial Computer Software - Use governed by
# terms of Orbital Insight commercial license agreement.

"""
Created on Tue Apr 17 14:37:23 2018

@author: fergal


"""
from __future__ import print_function
from __future__ import division

from pdb import set_trace as debug
import shutil
import sys
import os
import re

def main():
    #tmpdir is where I do my setup.
    #build_dir is where conda does its work. These must be different.
    #build_dir must be on a device that supports 256 character paths
    #(which encrytpted filesystems may not support)
    tmpdir = "./tmp/"  #<-- notice the dot
    build_dir = "./build"

    for path in [tmpdir, build_dir]:
        if os.path.exists(path):
            shutil.rmtree(path)
        os.mkdir(path)

    shutil.copy("./conda.yaml", "%s/conda.yaml" %(tmpdir))
    shutil.copy("./build.sh", "%s/build.sh" %(tmpdir))
#    shutil.copy("./setup.py", "%s/setup.py" %(tmpdir))

    try:
        conda_exe = os.environ['CONDA_EXE']
    except KeyError:
        print("CONDA_EXE env var not set. Don't know where conda binary is installed")
        return 1


    print("Running conda build...")
    cmd = "%s build %s --croot %s --output-folder ./channel" %(conda_exe, tmpdir, build_dir)
    status = os.system(cmd) #
    if status > 0:
        raise ValueError("conda build reported an error: Status: %i" %(status))


def get_git_package(giturl, package, tag, tmpdir):
    cwd = os.getcwd()
    os.chdir(tmpdir)
    status = os.system("git clone %s %s" %(giturl, package))

    if status > 0:
        os.chdir(cwd)
        raise ValueError("Git clone reported an error: Status: %i" %(status))

    os.chdir(package)
    status =  os.system("git checkout %s" %(tag))
    if status > 0:
        os.chdir(cwd)
        raise ValueError("Git checkout reported an error: Status: %i" %(status))

    os.chdir(cwd)

if __name__ == "__main__":
    main()
