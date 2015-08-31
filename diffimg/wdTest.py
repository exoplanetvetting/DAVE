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
import diffimg
import tpf


"""Test K2 difference imager on  a constant star
"""


def main():
    ar = mastio2.K2Archive()

    kepid = 206011438

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

    fits, hdr = ar.getTargetPixelFile(kepid, 3, header=True)
    cube = tpf.getTargetPixelArrayFromFits(fits, hdr)
    gain = hdr['gain']
    cube *= gain

    i0 = np.random.randint(0, len(cin))
#    i0 = 2722
    indexInTransit = [i0]

    mp.figure(1)
    try:
        diff, oot = diffimg.constructK2DifferenceImage(cube, indexInTransit, rollPhase, flags)
#        diffimg.plotDiffDiagnostic(diff, np.sqrt(np.fabs(oot)), np.sqrt(np.fabs(cube[i0])))
        diffimg.plotDiffDiagnostic(diff, oot, np.sqrt(np.fabs(oot)))
        mp.title('Poisson Noise')

    except ValueError, e:
        print e

#    mp.figure(2)
#    mp.clf()
#    rp0 = rollPhase[i0]
#    mp.axvline(time[i], color='r')
#    mp.axhline(rp0, color='grey')
#    mp.axvline(time[indexBefore], color='g')
#    mp.axvline(time[indexAfter], color='g')
#    mp.plot(time[~badIdx], rollPhase[~badIdx], 'k.')
#    mp.plot(time[~badIdx], rot[~badIdx, 1], 'r.')
#    plotThrusterFirings(flags,time)
