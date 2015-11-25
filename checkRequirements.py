# -*- coding: utf-8 -*-
"""
Created on Thu Nov 19 16:09:09 2015

@author: fergal
"""
from __future__ import print_function

import matplotlib as mp
import numpy as np

import subprocess
import pip


def checkInstalled(cmd):

    print("Checking for %s" %(cmd))
    try:
        res = subprocess.check_output(cmd.split())
    except OSError:
        print("... not found")
        return False

    print("... OK")
    return True


def checkImport(package, installIfNeeded=False):

    print("Checking for python package %s" %(package))
    try:
        __import__(package)
        print("... OK")
    except ImportError:
        print("... not found")

        if installIfNeeded:
            print("Attempting to install with pip...")
            val = pip.main(['install', package])
            if val > 0:
                print("Failed")
                return False
            else:
                print("Success")
                return True
        else:
            return False
    return True

def main():

    isOk = True
    isOk &= checkInstalled("python2.7 --version")
    isOk &= checkInstalled("g++ -dumpversion")
    isOk &= checkInstalled("octave --version")
#    isOk &= checkInstalled("")

    isOk &= checkImport("numpy")
    isOk &= checkImport("scipy")
    isOk &= checkImport("matplotlib")

    #This is specialised code.
    isOk &= checkImport("bls", True)

    print("*****************\n")
    if isOk:
        print("All requirements met!")
    else:
        print("Sorry, not all requirements met. Please check the output to see")
        print("which software you must install before installing Dave.")


if __name__ == "__main__":
    main()
