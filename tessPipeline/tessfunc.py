"""
Created on Tue Nov 27 21:44:04 2018

@author: fergal
"""

from __future__ import division

from pdb import set_trace as debug
import matplotlib.pyplot as plt
import matplotlib.patches as mpatch
import matplotlib as mpl
import pandas as pd
import numpy as np
from astroquery.mast import Tesscut
from astropy.coordinates import SkyCoord
from dave.tessPipeline.tessio import TessDvtLocalArchive
from dave.tessPipeline.tessmastio import TessAstroqueryArchive
import os

def serve(sector, tic, planetNum, localPath, source_):
    
    if planetNum > 1:
        planetNum = 1

    ar = TessAstroqueryArchive(localPath)
#    ar = TessDvtLocalArchive(localPath)
        
    if source_ == "tess_2min":
        dvt, hdr = ar.getDvt(tic, sector, ext=planetNum, header=True)
#    dvt, hdr = ar.getLightcurve(tic, sector, ext=planetNum, header=True)
        tpf, hdr_tpf = ar.getTPF(tic, sector, ext=planetNum, header=True)

    elif source_ == "tess_FFI":
        aperture_size = 11# in pixels
        tpf, hdr_tpf = ar.getFfiCutout(tic, sector, aperture_size)
        dvt, hdr = tpf, hdr_tpf

    return dvt, hdr, tpf, hdr_tpf

def serve_TESS_FFI(sector, tic, planetNum, localPath):

    ar = TessAstroqueryArchive(localPath)
    tpf, hdr_tpf = ar.getFfiCutout(tic, sector, 11)
    dvt, hdr = tpf, hdr_tpf
#    print(tpf_header)
#    print(tpf.shape)
#    xxxx

    return dvt, hdr, tpf, hdr_tpf

def getOutputBasename(basePath, tic):
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

    ticStr = str(int(tic))
    path = os.path.join(basePath, ticStr)

    if not os.path.exists(path):
        os.mkdir(path)
    if not (os.access(path, os.W_OK)):
        raise IOError("Can't write to output directory %s" %path)

    return os.path.join(path, ticStr)
