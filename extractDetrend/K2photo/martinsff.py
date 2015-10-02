import sys, os
import scipy
from pylab import *
from matplotlib import *
from scipy.stats import *
from numpy import *
from scipy import *
import kepfit
import kepmsg


"""
This code is based on the PyKE routine kepsff
found at keplerscience.arc.nasa.gov

The kepsff code is based on Vanderberg and Johnson 2014.
If you use this you must cite V&J 2014.

"""

def martinsff(intime,indata,centr1,centr2,
    npoly_cxcy,sigma_cxcy,npoly_ardx,
    npoly_dsdt,sigma_dsdt,npoly_arfl,sigma_arfl,verbose,logfile,
    status):

# startup parameters

    status = 0
    labelsize = 16
    ticksize = 14
    xsize = 20
    ysize = 8
    lcolor = '#0000ff'
    lwidth = 1.0
    fcolor = '#ffff00'
    falpha = 0.2
    seterr(all="ignore")





# fit centroid data with low-order polynomial

    cfit = zeros((len(centr2)))
    csig = zeros((len(centr2)))
    functype = 'poly' + str(npoly_cxcy)
    pinit = array([nanmean(centr2)])
    if npoly_cxcy > 0:
        for j in range(npoly_cxcy):
            pinit = append(pinit,0.0)
    try:
        coeffs, errors, covar, iiter, sigma, chi2, dof, fit, plotx, ploty, status = \
            kepfit.lsqclip(functype,pinit,centr1,centr2,None,sigma_cxcy,sigma_cxcy,10,logfile,verbose)
        for j in range(len(coeffs)):
            cfit += coeffs[j] * numpy.power(centr1,j)
            csig[:] = sigma
    except:
        message  = 'ERROR -- KEPSFF: could not fit centroid data with polynomial. There are no data points within the range of input rows %d - %d. Either increase the stepsize (with an appreciation of the effects on light curve quality this will have!), or better yet - cut the timeseries up to remove large gaps in the input light curve using kepclip.' % (t1,t2)
        status = kepmsg.err(logfile,message,verbose)
#        sys.exit('')
        os._exit(1)

# reject outliers

    time_good = array([],'float64')
    centr1_good = array([],'float32')
    centr2_good = array([],'float32')
    flux_good = array([],'float32')
    cad_good = array([],'int')
    for i in range(len(cfit)):
        if abs(centr2[i] - cfit[i]) < sigma_cxcy * csig[i]:
            time_good = append(time_good,intime[i])
            centr1_good = append(centr1_good,centr1[i])
            centr2_good = append(centr2_good,centr2[i])
            flux_good = append(flux_good,indata[i])

# covariance matrix for centroid time series

    centr = concatenate([[centr1_good] - mean(centr1_good), [centr2_good] - mean(centr2_good)])
    covar = cov(centr)

# eigenvector eigenvalues of covariance matrix

    [eval, evec] = numpy.linalg.eigh(covar)
    ex = arange(-10.0,10.0,0.1)
    epar = evec[1,1] / evec[0,1] * ex
    enor = evec[1,0] / evec[0,0] * ex
    ex = ex + mean(centr1)
    epar = epar + mean(centr2_good)
    enor = enor + mean(centr2_good)

# rotate centroid data

    centr_rot = dot(evec.T,centr)

# fit polynomial to rotated centroids

    rfit = zeros((len(centr2)))
    rsig = zeros((len(centr2)))
    functype = 'poly' + str(npoly_ardx)
    pinit = array([nanmean(centr_rot[0,:])])
    pinit = array([1.0])
    if npoly_ardx > 0:
        for j in range(npoly_ardx):
            pinit = append(pinit,0.0)
    try:
        coeffs, errors, covar, iiter, sigma, chi2, dof, fit, plotx, ploty, status = \
            kepfit.lsqclip(functype,pinit,centr_rot[1,:],centr_rot[0,:],None,100.0,100.0,1,
                           logfile,verbose)
    except:
        message  = 'ERROR -- KEPSFF: could not fit rotated centroid data with polynomial'
        status = kepmsg.err(logfile,message,verbose)
    rx = linspace(nanmin(centr_rot[1,:]),nanmax(centr_rot[1,:]),100)
    ry = zeros((len(rx)))
    for i in range(len(coeffs)):
        ry = ry + coeffs[i] * numpy.power(rx,i)

# calculate arclength of centroids

    s = zeros((len(rx)))
    for i in range(1,len(s)):
        work3 = ((ry[i] - ry[i-1]) / (rx[i] - rx[i-1]))**2
        s[i] = s[i-1] + math.sqrt(1.0 + work3) * (rx[i] - rx[i-1])

# fit arclength as a function of strongest eigenvector

    sfit = zeros((len(centr2)))
    ssig = zeros((len(centr2)))
    functype = 'poly' + str(npoly_ardx)
    pinit = array([nanmean(s)])
    if npoly_ardx > 0:
        for j in range(npoly_ardx):
            pinit = append(pinit,0.0)
    try:
        acoeffs, errors, covar, iiter, sigma, chi2, dof, fit, plotx, ploty, status = \
            kepfit.lsqclip(functype,pinit,rx,s,None,100.0,100.0,100,logfile,verbose)
    except:
        message  = 'ERROR -- KEPSFF: could not fit rotated centroid data with polynomial'
        status = kepmsg.err(logfile,message,verbose)

# correlate arclength with detrended flux

    t = copy(time_good)
    y = copy(flux_good)
    z = centr_rot[1,:]
    x = zeros((len(z)))
    for i in range(len(acoeffs)):
        x = x + acoeffs[i] * numpy.power(z,i)

# calculate time derivative of arclength s

    dx = zeros((len(x)))
    for i in range(1,len(x)):
        dx[i] = (x[i] - x[i-1]) / (t[i] - t[i-1])
    dx[0] = dx[1]

# fit polynomial to derivative and flag outliers (thruster firings)

    dfit = zeros((len(dx)))
    dsig = zeros((len(dx)))
    functype = 'poly' + str(npoly_dsdt)
    pinit = array([nanmean(dx)])
    if npoly_dsdt > 0:
        for j in range(npoly_dsdt):
            pinit = append(pinit,0.0)
    try:
        dcoeffs, errors, covar, iiter, dsigma, chi2, dof, fit, dumx, dumy, status = \
            kepfit.lsqclip(functype,pinit,t,dx,None,3.0,3.0,10,logfile,verbose)
    except:
        message  = 'ERROR -- KEPSFF: could not fit rotated centroid data with polynomial'
        status = kepmsg.err(logfile,message,verbose)
    for i in range(len(dcoeffs)):
        dfit = dfit + dcoeffs[i] * numpy.power(t,i)
    centr1_pnt = array([],'float32')
    centr2_pnt = array([],'float32')
    time_pnt = array([],'float64')
    flux_pnt = array([],'float32')
    dx_pnt = array([],'float32')
    s_pnt = array([],'float32')
    time_thr = array([],'float64')
    flux_thr = array([],'float32')
    dx_thr = array([],'float32')
    thr_cadence = zeros(len(t),dtype=bool)
    for i in range(len(t)):
        if dx[i] < dfit[i] + sigma_dsdt * dsigma and dx[i] > dfit[i] - sigma_dsdt * dsigma:
            time_pnt = append(time_pnt,time_good[i])
            flux_pnt = append(flux_pnt,flux_good[i])
            dx_pnt = append(dx_pnt,dx[i])
            s_pnt = append(s_pnt,x[i])
            centr1_pnt = append(centr1_pnt,centr1_good[i])
            centr2_pnt = append(centr2_pnt,centr2_good[i])
        else:
            time_thr = append(time_thr,time_good[i])
            flux_thr = append(flux_thr,flux_good[i])
            dx_thr = append(dx_thr,dx[i])
            thr_cadence[i] = True

# fit arclength-flux correlation

    cfit = zeros((len(time_pnt)))
    csig = zeros((len(time_pnt)))
    functype = 'poly' + str(npoly_arfl)
    pinit = array([nanmean(flux_pnt)])
    if npoly_arfl > 0:
        for j in range(npoly_arfl):
            pinit = append(pinit,0.0)
    try:
        ccoeffs, errors, covar, iiter, sigma, chi2, dof, fit, plx, ply, status = \
            kepfit.lsqclip(functype,pinit,s_pnt,flux_pnt,None,sigma_arfl,sigma_arfl,100,logfile,verbose)
    except:
        message  = 'ERROR -- KEPSFF: could not fit rotated centroid data with polynomial'
        status = kepmsg.err(logfile,message,verbose)

# correction factors for unfiltered data

    centr = concatenate([[centr1] - mean(centr1_good), [centr2] - mean(centr2_good)])
    centr_rot = dot(evec.T,centr)
    yy = copy(indata)
    zz = centr_rot[1,:]
    xx = zeros((len(zz)))
    cfac = zeros((len(zz)))
    for i in range(len(acoeffs)):
        xx = xx + acoeffs[i] * numpy.power(zz,i)
    for i in range(len(ccoeffs)):
        cfac = cfac + ccoeffs[i] * numpy.power(xx,i)

# apply correction to flux time-series

    out_detsap = indata / cfac

    return out_detsap, cfac, thr_cadence

