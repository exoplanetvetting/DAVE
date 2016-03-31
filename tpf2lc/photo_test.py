from __future__ import division, print_function

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits as pyfits


## from photutils.aperture import aperture_elliptical

def raw_moment(data, iord, jord):
    nrows, ncols = data.shape
    y, x = np.mgrid[:nrows, :ncols]
    data = data * x**iord * y**jord
    return data.sum()

def intertial_axis(data):
    """Calculate the x-mean, y-mean, and cov matrix of an image."""
    data_sum = data.sum()
    m10 = raw_moment(data, 1, 0)
    m01 = raw_moment(data, 0, 1)
    x_bar = m10 / data_sum
    y_bar = m01 / data_sum
    u11 = (raw_moment(data, 1, 1) - x_bar * m01) / data_sum
    u20 = (raw_moment(data, 2, 0) - x_bar * m10) / data_sum
    u02 = (raw_moment(data, 0, 2) - y_bar * m01) / data_sum
    cov = np.array([[u20, u11], [u11, u02]])
    return x_bar, y_bar, cov

def plot_bars(x_bar, y_bar, cov, ax):
    """Plot bars with a length of 2 stddev along the principal axes."""
    def make_lines(eigvals, eigvecs, mean, i):
        """Make lines a length of 2 stddev."""
        std = np.sqrt(eigvals[i])
        vec = 3 * std * eigvecs[:,i] / np.hypot(*eigvecs[:,i])
        x, y = np.vstack((mean-vec, mean, mean+vec)).T
        return x, y
    mean = np.array([x_bar, y_bar])
    eigvals, eigvecs = np.linalg.eigh(cov)
    ax.plot(*make_lines(eigvals, eigvecs, mean, 0), marker='o', color='white')
    ax.plot(*make_lines(eigvals, eigvecs, mean, -1), marker='o', color='red')
    ax.axis('image')

def polyfit_iter(time,flux,order):
    pass


"""
if __name__ == '__main__':
    dirname = '/Users/tom/Projects/K2_science/Jan2014_preSafeMode/LC'
    filename = '/kplr060020058-2014020013418_lpd-targ.fits'
    f = pyfits.open(dirname + filename)
    time = f[1].data['TIME'] - f[1].data['TIME'][0]
    fluxarr = f[1].data['FLUX']

    tmask = np.isfinite(time)
    time, fluxarr = time[tmask], fluxarr[tmask,:,:]

    #pixmaskx = (18,33)
    #pixmasky = (17,33)

    pixmaskx = (0,50)
    pixmasky = (0,50)

    fluxbox = fluxarr[:,pixmasky[0]:pixmasky[1],pixmaskx[0]:pixmaskx[1]]

    flatim = np.sum(fluxbox,axis=0)

    brightestpix = np.unravel_index(flatim.argmax(), flatim.shape)
    bpixy, bpixx = brightestpix
    bpixy = np.where(bpixy<0,0,bpixy)
    bpixx = np.where(bpixx<0,0,bpixx)

    win = 17
    pixmaskx = (bpixx-win,bpixx+win+1)
    pixmasky = (bpixy-win,bpixy+win+1)

    fluxbox = fluxarr[:,pixmasky[0]:pixmasky[1],pixmaskx[0]:pixmaskx[1]]
    flatim = np.sum(fluxbox,axis=0)

    flatnormim = (flatim / np.median(flatim)) - 1.0

    #xbar and ybar are the means
    xbar, ybar, cov = intertial_axis(flatnormim)
    mean = np.array([xbar, ybar])
    eigvals, eigvecs = np.linalg.eigh(cov)

    do_plot = True
    if do_plot:
        fig, ax = plt.subplots()
        ax.imshow(flatnormim,interpolation='nearest',
            origin='lower')
        plot_bars(xbar, ybar, cov, ax)

    nsigma=12.
    stdevb_fix = np.sqrt(eigvals[-1])
    vec = nsigma * stdevb_fix * eigvecs[:,-1] / np.hypot(*eigvecs[:,-1])
    angle_fix = np.arctan(vec[0] / vec[1]) - (0.5*np.pi)

    stdeva_fix = np.sqrt(eigvals[0])

    ph = aperture_elliptical(flatnormim,xbar,ybar,nsigma * stdevb_fix,
        nsigma * stdeva_fix,angle_fix,method='subpixel',)


    xbar_fix, ybar_fix, cov_fix = xbar, ybar, cov
    lc = np.zeros_like(time)
    xpix = np.zeros_like(time)
    ypix = np.zeros_like(time)
    for i in range(np.shape(fluxbox)[0]):
        if np.all(fluxbox[i,:,:] == 0):
            continue
        normim = (fluxbox[i,:,:] / np.median(fluxbox[i,:,:])) - 1.0
        xbar, ybar, cov = intertial_axis(normim)
        mean = np.array([xbar, ybar])
        eigvals, eigvecs = np.linalg.eigh(cov)
        stdevb = np.sqrt(eigvals[-1])
        vec = nsigma * stdevb * eigvecs[:,-1] / np.hypot(*eigvecs[:,-1])
        angle = np.arctan(vec[0] / vec[1]) - (0.5*np.pi)

        stdeva = np.sqrt(eigvals[0])


        ph = aperture_elliptical(normim,xbar,ybar,nsigma * stdevb_fix,
            nsigma * stdeva_fix,angle_fix,method='exact',subpixels=25)
        lc[i] = ph
        xpix[i] = xbar
        ypix[i] = ybar

    ## lets fit a quadratic function
    time,lc = time[lc > 0], lc[lc>0]
    lc  = (lc /np.median(lc)) - 1.0


    #i = 0
    #m2 = np.ones_like(time,dtype=bool)
    #while i < 5:
    #    polyf = np.polyfit(time[m2], lc[m2], deg=3)
    #    tfit = np.polyval(polyf,time[m2])
    #    resid = lc[m2] - tfit
    #    mask = np.abs(resid) <= 5.* np.std(resid)
    #    m2 = mask
    #    i += 1

    from scipy.ndimage.filters import median_filter
    window = 48
    medval = median_filter(lc,window,mode='mirror')

    #sigmaclip
    mlc = lc - medval
    sig = np.std(mlc)
    medm = (mlc < 5*sig)& (mlc > -5*sig)
    t1,f1 = time[medm], mlc[medm]
    #from scipy.ndimage import generic_filter
    time,f_norm = t1,f1
    freqs = np.linspace((1. / 1.1) * (2.* np.pi), (30*24.) * (2.*np.pi),50000)
    from scipy.signal import lombscargle as lomb
    l = lomb(time,f_norm,freqs)
    fig = plt.figure(figsize=[9,7])
    ax1 = fig.add_subplot(211)
    ax1.plot(time,f_norm*1.E2,'k.-')
    ax1.set_xlabel('Time (days)',fontsize=16)
    ax1.set_ylabel('A (\%)',fontsize=16)
    ax1.set_xlim([min(time),max(time)])
    ax1.set_ylim([-7,7])
    ax2 = fig.add_subplot(212)
    ax2.plot(freqs / (2.*np.pi) / 86400. * 1.E6,
                 np.sqrt(4*(l/time.shape[0])) * 1.E2,'k-')
    ax2.set_xlabel('\mu Hz',fontsize=16)
    ax2.set_ylabel('A (\%)',fontsize=16)
    ax2.set_xlim([450,2500])
    ax2.set_ylim([0,0.9])
    #ax1.set_title('GD1212',fontsize=18)
    plt.minorticks_on()
    plt.tight_layout()






"""






