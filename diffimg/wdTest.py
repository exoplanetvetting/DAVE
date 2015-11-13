# -*- coding: utf-8 -*-
"""
Created on Sat Aug 29 21:07:11 2015

@author: fergal

$Id$
$URL$
"""

import matplotlib.pyplot as mp
import numpy as np

import mastio2
import play as roll
import diffimg2 as diffimg
import tpf


"""Test K2 difference imager on  a constant star
"""


def main():
    ar = mastio2.K2Archive()

    kepid = 206011438 #WD
    kepid = 206103150 #WASP star

    fits = ar.getLongCadence(kepid, 3)
    time = fits['TIME']
    cin = fits['CADENCENO']
    flags = fits['SAP_QUALITY']
    pa = fits['SAP_FLUX']
    pdc = fits['PDCSAP_FLUX']
    cent1 = fits['MOM_CENTR1']
    cent2 = fits['MOM_CENTR2']
    badIdx = flags > 0

    #Compute roll phase
    centColRow = np.vstack((cent1, cent2)).transpose()
    rot = roll.computeArcLength(centColRow, flags>0)
    rollPhase = rot[:,0]
    rollPhase[flags>0] = -9999    #A bad value

    fits, hdr = ar.getLongTpf(kepid, 3, header=True)
    cube = tpf.getTargetPixelArrayFromFits(fits, hdr)
    gain = hdr['gain']
    cube *= gain

#    i0 = np.random.randint(0, len(cin))
    i0 = 2622

    mp.figure(1)
    itr = cube[i0]
    oot = diffimg.getInterpolatedOotImage(cube, rollPhase, flags, i0)
    diff = itr-oot
    diffimg.plotDiffDiagnostic(diff, oot, itr)
