
#Tools for plotting target pixel files
#$Id: plotTpf.py 2159 2015-09-22 18:26:13Z fmullall $
#$URL: svn+ssh://flux/home/fmullall/svn/kepler/py/plotTpf.py $
#Fergal Mullally


from  matplotlib.widgets import MultiCursor
import matplotlib.pyplot as mp
import numpy as np
import dave.fileio.tpf as tpf


def plotCadence(singleImage, hdr=None, axis='Fixed', *args, **kwargs):
    """
    Plot a single cadence from a TPF data cube

    Inputs:
    ---------
    singleImage
        (np 2d array) Image to plot

    Optional Inputs:
    ----------------
    hdr
        Header from fits file
    axis
        ('relative' of 'fixed'). If set to fixed, and hdr is set,
        the axis labels are in units of pixels. If relative,
        the bottom left of the image is labeled as pixel 0,0

    Additional arguments passed to mp.imshow
    """
    #Default interpolation is to plot square pixels
    interpolation = kwargs.pop('interpolation', 'nearest')
    cmap = kwargs.pop('cmap', mp.cm.YlGnBu_r)

    shape = singleImage.shape   #Gives the shape as (row,col)

    #What numbers to put on x and y axes?
    extent=None
    if axis.lower()=="relative" or hdr is None:
        nr, nc = singleImage.shape
        extent=[0,nc, 0, nr]
    elif axis.lower() == "fixed":

        #For kepler TPF files, all the columns in the binary
        #table contain images of the same shape, so it
        #doesn't matter which CRPIXs I take.
        c0 = float(hdr['1CRV4P'])
        r0 = float(hdr['2CRV4P'])

        extent = [c0, c0+shape[1], r0, r0+shape[0]]
    else:
        raise ValueError("Unrecognised axis style requested")

    mp.imshow(singleImage, origin='lower', interpolation=interpolation, \
        extent=extent, cmap=cmap, *args, **kwargs)



def plotOptimalAperture(mask, *args, **kwargs):
    """Won't work if img is plotted  in ra/dec space
    """
    nR, nC = mask.shape
    ax = mp.gca()
    x0, x1, y0, y1 = mp.axis()

    for r in range(nR):
        for c in range(nC):
            if mask[r, c] == 3:
                sq = mp.Rectangle([c+x0, r+y0], 1,1, \
                    edgecolor='r', lw=2, fill=False)
                ax.add_patch(sq)




def plotTpfCadenceDiagnostic(index, fits, hdr):
    mp.clf()

    mp.subplot(221)
    cube = tpf.getTargetPixelArrayFromFits(fits, hdr, 'RAW_CNTS')
    #Missing values in raw counts are set to -1.
    #We want them set to NaN, like everything else
    idx = cube[index, :,:] < 0
    cube[index, idx] = np.nan

    plotCadence(cube[index,:,:], hdr)
    mp.title('Raw Counts')
    mp.colorbar()

    mp.subplot(222)
    cube = tpf.getTargetPixelArrayFromFits(fits, hdr, 'FLUX')
    plotCadence(cube[index,:,:], hdr)
    mp.title('Flux (e/s)')
    mp.colorbar()

    mp.subplot(223)
    cube = tpf.getTargetPixelArrayFromFits(fits, hdr, 'COSMIC_RAYS')
    plotCadence(cube[index,:,:], hdr)
    mp.title('Cosmic Rays (e/s)')
    mp.colorbar()

    mp.subplot(224)
    cube = tpf.getTargetPixelArrayFromFits(fits, hdr, 'FLUX_BKG')
    plotCadence(cube[index,:,:], hdr)
    mp.title('Background (e/s)')
    mp.colorbar()


    #kepid = hdr['keplerid']
    #quarter = 99 #hdr['quarter']
    #rin = index+1
    #cin = rin + fits['CADENCENO'][0]
    #title = "Q%i %i %s. CIN(RIN) = %i(%i)" % \
        #(quarter, kepid, "??", cin, rin)
    #mp.suptitle(title)


def plotTpfLc(cube, hdr, *args, **kwargs):
    """Plot a pixel time series for every pixel in a cube

    Similar, static plots, used to be included in the DV report
    but were removed for some reason.

    Inputs:
    --------
    cube    (np 3d array)
        A TPF data cube as produced by tpf.getTargetPixelArrayFromFits
    hdr     (FITS header)
        The header of the TPF file


    Optional Inputs:
    --------------
    detrend     (boolean)
        Whether to apply a Savitzy-Golay filter to the data
        Note, transits are not protected, so this will distort
        your lightcurve in funny ways. Default: False

    big (boolean)
        Force plotting of really big masks (greater than 49 pixels)
        Default: False

    All other optional arguments are passed to mp.plot()

    Returns:
    ----------
    **None**

    Todo:
    --------
    * Optional argument to set the time axis
    * Show axis ticks on the bottom and left plots
    * Ability to highlight certain time intervals (e.g times
      of transits
    * Ability to specify detrending algorithm.

    """

    detrendFlag = kwargs.pop('detrend', False)
    bigFlag = kwargs.pop('big', False)

    nCin, nR, nC = cube.shape

    print(nCin, nR, nC)
    nPix = nR*nC

    #Prevent automatic plotting of really big masks.
    if nPix > 49 and not bigFlag:
        raise ValueError("Too many pixels. Set big=True to force plotting")

    ax0 = mp.subplot(nR, nC, 1)
    mp.gcf().subplots_adjust(wspace=0, hspace=0)

    mny = np.inf
    mxy = -np.inf
    time = np.arange(nCin)
    for i in range(nPix):
        ax = mp.subplot(nR, nC, i+1, sharex=ax0, sharey=ax0)
        ax.set_xticks([])
        ax.set_yticks([])

        col = np.fmod(i, nC)
        row = nR - 1 - np.floor( (i)/float(nC))
#        print i+1, col, row

        if col==0:
            mp.ylabel("%i" % row, rotation="horizontal")

        if row==0:
            mp.xlabel("%i" %col)
        lc = cube[:, row, col]
        if np.any(np.isfinite(lc)) and 1:
            idx= np.isfinite(lc)
            lc /= np.median(lc[idx])
            lc -= 1

            if detrendFlag:
                gapFilledFlux = sgolay.fillGaps(time, lc, ~idx)
                lowPass = ss.savgol_filter(gapFilledFlux, 48+1, 3, mode='nearest')
                lc -= lowPass

            mp.plot(time[idx], lc[idx], *args, **kwargs)

            mny = min(mny, np.min(lc[idx]))
            mxy = max(mxy, np.max(lc[idx]))

    mp.ylim(mny, mxy)
