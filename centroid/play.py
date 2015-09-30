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
    quarter = 14

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
    col += bbox[0] + .5
    row += bbox[2] + .5

#    col,row = 670.5, 653.5

    print "Initial Guess for Q%i (%.3f %.3f)" %(quarter, col, row)
#    x0 = [120e3]
#    args = (module, output, col, row, bbox, img)
#    bounds=[(0, None)]
#    res = sopt.minimize(costFunc1, x0, args, method="L-BFGS-B", \
#        bounds=bounds)

    x0 = [col, row, 120e3]
    scale = 1 #np.max(img)
    args = (module, output, bbox, img/scale)
    options = {'disp':False, 'eps':.02, 'maxiter':80}
    bounds=[(bbox[0], bbox[1]), (bbox[2], bbox[3]), (1, None)]
    res = sopt.minimize(costFunc2, x0, args, method="L-BFGS-B", \
        bounds=bounds, options=options)
    res.x[2] *= scale  #Rescale amplitude
#    res.x[:2] += .02

#    c0 = np.floor(col)
#    r0 = np.floor(row)
#    score = np.zeros((50,50))
#    for i in range(0, 50, 1):
#        c = c0 + i/50.
#        print i
#        for j in range(0, 50, 1):
#            r = r0 + j/50.
#            x0 = [120e3]
#            args = (module, output, c, r, bbox, img)
#            bounds=[(0, None)]
#            res = sopt.minimize(costFunc1, x0, args, \
#                method="L-BFGS-B", bounds=bounds)
#            score[j, i] = res.fun
#
#    mp.figure(2)
#    mp.clf()
#    mp.imshow(np.log10(score), cmap=mp.cm.rainbow, interpolation="nearest", origin="bottom")
#    mp.colorbar()
#    r,c = np.unravel_index(np.argmax(img), img.shape)
#    print c0+c, r0+r
#
#    return score


    mp.figure(1)
    mp.clf()
    mp.subplot(141)
    plotTpf.plotCadence(img, hdr)
    mp.colorbar()

    mp.subplot(142)
    kPrf = prf.KeplerPrf("/home/fergal/data/keplerprf")
    c,r = res.x[0], res.x[1]
    model = kPrf.getPrfForBbox(module, output, c, r, bbox)
    model *= res.x[2]
    plotTpf.plotCadence(model, hdr)
    mp.colorbar()

    mp.subplot(143)
    diff = img-model
    plotTpf.plotCadence(diff, hdr)
    mp.colorbar()

    print "Performance %.3f" %(np.max(np.abs(diff))/np.max(img))

    return res

def costFunc2(x, module, output, bbox, img):
    kPrf = prf.KeplerPrf("/home/fergal/data/keplerprf")
    model = kPrf.getPrfForBbox(module, output, x[0], x[1], bbox)
    model *= x[2]

    cost = img-model
    cost = np.sum(cost**2)
#    print "%.3e" %(cost)
#    cost = np.log10(cost)
    return cost


def costFunc1(x, module, output, col, row, bbox, img):
    kPrf = prf.KeplerPrf("/home/fergal/data/keplerprf")
    model = kPrf.getPrfForBbox(module, output, col, row, bbox)
    model *= x[0]

    cost = img-model
    cost = np.sum(cost**2)
    return cost



def main2():
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
#    func2 = lambda x: np.max(np.abs(func(x)))

    mp.figure(1)
    mp.clf()
    mp.subplot(141)
    plotTpf.plotCadence(img, hdr)
    mp.plot(col, row, 'ro')
    mp.colorbar()

    kPrf = prf.KeplerPrf("/home/fergal/data/keplerprf")
    model = kPrf.getPrfForBbox(module, output, col, row, bbox)

#    mp.subplot(142)
#    plotTpf.plotCadence(model, hdr, cmap=mp.cm.jet)
#    mp.colorbar()
#
#    mp.subplot(143)
    scale = np.max(img)/np.max(model)
#    diff = img - scale*model
#    plotTpf.plotCadence(diff, hdr)
#    print "M1", np.log10(np.sum(diff**2)), scale
#    print "M2", np.log10(func2([col, row, scale])), 3*np.max(img)
#    mp.colorbar()

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
    x0 = [col+.5, row+.5, scale]
#    x0 = [670.78, 653.87, scale]
    print x0, bbox
    method = "L-BFGS-B"
#    method = "SLSQP"
    result = sopt.minimize(func2, x0, method=method, bounds=[(bbox[0], bbox[1]), (bbox[2], bbox[3]), (0, None)]  )
    c, r = result.x[:2]
    c, r = 670.875, 653.98

    model = kPrf.getPrfForBbox(module, output, c, r, bbox)
    mp.subplot(142)
    mp.cla()
    plotTpf.plotCadence(model, hdr)
    mp.colorbar()

    mp.subplot(143)
    sModel = model * result.x[2]  #Scale model to flux
    plotTpf.plotCadence(sModel, hdr)
    mp.colorbar()

    mp.subplot(144)
    diff = img - sModel
    plotTpf.plotCadence(diff, hdr)
    mp.colorbar()

    return result


#    #globaal search
#    score = np.zeros((50,50))
#    c0, r0 = np.floor(col), np.floor(row)
#    for i in range(0, 50, 1):
#        print i
#        c = c0 + i/50.
#        for j in range(0, 50, 1):
#            r = r0 + j/50.
#
#            score[i, j] = np.log10(func2([c, r, scale]))
#
#    mp.figure(2)
#    mp.clf()
#    mp.subplot(121)
#    mp.imshow(score, origin="bottom", cmap=mp.cm.rainbow, interpolation="nearest")
#    mp.colorbar()
#
#    mp.subplot(122)
#    loc = np.unravel_index(np.argmin(score), (50,50))
#    c = np.floor(col) + loc[0]/50.
#    r = np.floor(row) + loc[1]/50.
#    model = kPrf.getPrfForBbox(module, output, c, r, bbox)
#    scale = np.max(img)/np.max(model)
#    diff = img - scale*model
#    plotTpf.plotCadence(diff, hdr)
#    mp.colorbar()
#
#    print c, r
#    print score[loc]
#    return score


def modelFunc(module, output, col, row, scale, bbox):
    kPrf = prf.KeplerPrf("/home/fergal/data/keplerprf")
    model = kPrf.getPrfForBbox(module, output, col, row, bbox)

    assert(np.all( np.isfinite(model)))
    vals = model.reshape((-1,1)).squeeze()
    return vals*scale