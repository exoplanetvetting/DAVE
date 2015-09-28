# -*- coding: utf-8 -*-
"""
Created on Fri Sep 25 08:59:40 2015

@author: fergal

$Id$
$URL$
"""

__version__ = "$Id$"
__URL__ = "$URL$"



import matplotlib.pyplot as mp
import numpy as np

import dave.fileio.mastio as mastio
import dave.fileio.tpf as tpf
import dave.centroid.prf as prf
import plotTpf

import scipy.optimize as sopt


def getExtent(img, hdr):
    shape= img.shape
    c0 = float(hdr['1CRV4P'])
    r0 = float(hdr['2CRV4P'])


    extent = [c0, c0+shape[1], r0, r0+shape[0]]
    return extent



def main():
    kepid = 8554498
    quarter = 4
    #I need to find these out blindly, but good enough for now
#    col, row = 670.74, 653.94  #Q4
#    col, row = 687.6, 639.74  #Q3
#    col, row = 684.58, 643.38  #Q1
#    col, row = 673.77 - .5 , 658.443 -.5


    ar = mastio.KeplerArchive()
    fits, hdr = ar.getLongTpf(kepid, quarter, header=True)
    hdr0 = ar.getLongTpf(kepid, quarter, ext=0)
    cube = tpf.getTargetPixelArrayFromFits(fits, hdr)

    module = hdr0['MODULE']
    output = hdr0['OUTPUT']
    img = cube[100]
    idx = np.isfinite(img)
    img[~idx] = 0
    bbox = getExtent(img, hdr)
    row, col = np.unravel_index(np.argmax(img), img.shape)
    col += bbox[0]
    row += bbox[2]

#    col, row = 687.6, 639.74
#    col, row = 687.82, 639.44
    print col, row
    print bbox

    #Score functions
    func = lambda x: img.reshape((-1,1)).squeeze() -\
        modelFunc(module, output, x[0], x[1], x[2], bbox)
    func2 = lambda x: np.sum(func(x)**2)


    mp.figure(1)
    mp.clf()
    mp.subplot(131)
    plotTpf.plotCadence(img, hdr, cmap=mp.cm.jet)
    mp.plot(col, row, 'ro')
    mp.colorbar()

    kPrf = prf.KeplerPrf("/home/fergal/data/keplerprf")
    model = kPrf.getPrfForBbox(module, output, col, row, bbox)

    mp.subplot(132)
    plotTpf.plotCadence(model, hdr, cmap=mp.cm.jet)
    mp.colorbar()

    mp.subplot(133)
    scale = np.max(img)/np.max(model)
    diff = img - scale*model
    plotTpf.plotCadence(diff, hdr)
    print "M1", np.log10(np.sum(diff**2)), scale
    print "M2", np.log10(func2([col, row, scale])), 3*np.max(img)
    mp.colorbar()

#    return
##    Plot a sequence of models while you vary row/col
#    for i in range(1,6):
#        mp.subplot(1, 5, i)
#        c = np.floor(col) + .5 +  .2*i
#        r = np.floor(row)
#        model = kPrf.getPrfAtColRow(module, output, c, r)
#        plotTpf.plotCadence(model, hdr)
#        mp.title("|C,r> = |%.2f %.2f>" %(c, r))
#    return


#    Fitting
    x0 = [col, row, scale]
#    result = sopt.leastsq(func, x0, full_output=True)
    result = sopt.minimize(func2, x0, method="L-BFGS-B", bounds=[(bbox[0], bbox[1]), (bbox[2], bbox[3]), (0, None)]  )
    return result


    #globaal search
    score = np.zeros((50,50))
    c0, r0 = np.floor(col), np.floor(row)
    for i in range(0, 50, 1):
        print i
        c = c0 + i/50.
        for j in range(0, 50, 1):
            r = r0 + j/50.

            score[i, j] = np.log10(func2([c, r, scale]))

    mp.figure(2)
    mp.clf()
    mp.subplot(121)
    mp.imshow(score, origin="bottom", cmap=mp.cm.rainbow, interpolation="nearest")
    mp.colorbar()

    mp.subplot(122)
    loc = np.unravel_index(np.argmin(score), (50,50))
    c = np.floor(col) + loc[0]/50.
    r = np.floor(row) + loc[1]/50.
    model = kPrf.getPrfForBbox(module, output, c, r, bbox)
    scale = np.max(img)/np.max(model)
    diff = img - scale*model
    plotTpf.plotCadence(diff, hdr)
    mp.colorbar()

    print c, r
    print score[loc]
    return score


def modelFunc(module, output, col, row, scale, bbox):
    kPrf = prf.KeplerPrf("/home/fergal/data/keplerprf")
    model = kPrf.getPrfForBbox(module, output, col, row, bbox)

    assert(np.all( np.isfinite(model)))
    vals = model.reshape((-1,1)).squeeze()
    return vals*scale