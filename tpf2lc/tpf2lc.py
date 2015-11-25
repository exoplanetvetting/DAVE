#Tools for manipulating target pixel files
#Tom Barclay
#Fergal Mullally
#Knicole Colon (edited 2015-10-29)


#Light Curve Extraction from TPFs (Optimal Apertures)
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import mastio

from astropy.io import fits as pyfits
from photo_test import raw_moment, intertial_axis, plot_bars

from astropy.stats.funcs import median_absolute_deviation as MAD
from scipy.ndimage import label
#from photo_test import raw_moment, intertial_axis, plot_bars


def getTargetPixelArrayFromFits(fits, hdr, column='FLUX'):
    nRows = int(hdr['NAXIS2'])
    #All imgs are the same size, so it doesn't matter which TDIM
    #I pick here.
    shape = eval(hdr['TDIM5'])
    #The reverse puts the columns on the x-axis and row
    #on the y-axis when you call mp.imshow()
    shape = tuple(reversed(shape))

    tpfArray = np.empty((nRows, shape[0], shape[1]))

    for i, cadence in enumerate(fits[column]):
        tpfArray[i,:,:] = cadence.reshape(shape)
        
    return tpfArray


def lightCurveFromTpfCube(cube, mask=None):
    """Extract a lightcurve from a TPF cube.

    Inputs:
    cube    (3d numpy array) As returned by getTargetPixelArrayFromFits()
    mask    (2d numpy array) Which pixels to use in extraction. Default
            is all finite pixels

    Returns
    A 1d numpy array

    Notes:
    Pixels marked Nan (i.e outside the aperture) are ignored regardless
    of the mask value.

    Default mask includes all finite pixels.

    To extract only the optimal aperture use
    mask = np.where(mask == 3, 1, 0)
    lightCurveFromTpfCube(cube, mask)
    """

    if mask is None:
        mask = np.ones( cube[0,:,:].shape)

    n, nR, nC = cube.shape
    out = np.zeros(n)
    for i in range(n):
        img = cube[i, :, :]
        idx = np.bitwise_and(mask > 0, np.isfinite(img))
        out[i] = np.sum(img[idx])

    return out


def bg_sub(fla):
    """
    subtract the background from a series of images
    by assuming the aperture is large enough to be
    predominantly background
    """
    for i in xrange(np.shape(fla)[0]):
        fla[i,:,:] = fla[i,:,:] - np.nanmedian(fla[i,:,:])
    return fla


def optimalAperture(t_time, t_fluxarr, t_quality, qual_cut=False, return_qual=False, toss_resat=True, bg_cut=5, skip_cads=None, skip=None):

    if skip is not None:
        skip = skip_cads
    else:
        skip = 0

    
    time = t_time[skip:]
    fluxarr = t_fluxarr[skip:]
    quality = t_quality[skip:]

    if qual_cut:
        time = time[quality == 0]
        fluxarr = fluxarr[quality == 0,:,:]
    elif toss_resat:
        # cadences where there is a wheel resaturation event
        time = time[quality != 32800]
        fluxarr = fluxarr[quality != 32800,:,:]

    #remove any nans
    fluxarr[fluxarr == 0] = np.nan

    #subtract background
    flux_b = bg_sub(fluxarr)

    # create a median image to calculate where the pixels to use are
    flatim = np.nanmedian(flux_b,axis=0)

    #find pixels that are X MAD above the median
    vals = flatim[np.isfinite(flatim)].flatten()
    mad_cut = 1.4826 * MAD(vals) * bg_cut

    region = np.where(flatim > mad_cut,1,0)
    lab = label(region)[0]

    #find the central pixel
    imshape = np.shape(flatim)
    centralpix = [1+imshape[0] // 2,1+imshape[1] // 2]

    #find brightest pix within 9x9 of central pix
    #this assumes target is at center of postage stamp which I think is ok
    centflatim = flatim[centralpix[0]-4:centralpix[0]+4,
                centralpix[1]-4:centralpix[1]+4]
    flatimfix = np.where(np.isfinite(centflatim),centflatim,0)
    brightestpix = np.unravel_index(flatimfix.argmax(), centflatim.shape)
    bpixy, bpixx = brightestpix


    #use all pixels in the postage stamp that are X MAD above the median
    #this identifies location of brightest pixel only
    regnum = lab[centralpix[0]-4+bpixy,centralpix[1]-4+bpixx]
    
    lc = np.zeros_like(time)
    xbar = np.zeros_like(time)
    ybar = np.zeros_like(time)

    #make a rectangular aperture for the moments thing
    ymin = np.min(np.where(lab == regnum)[0])
    ymax = np.max(np.where(lab == regnum)[0])
    xmin = np.min(np.where(lab == regnum)[1])
    xmax = np.max(np.where(lab == regnum)[1])
    #print(np.where(lab == regnum))
    momlims = [ymin,ymax+1,xmin,xmax+1]

    
    #loop that performs the aperture photometry
    for i,fl in enumerate(fluxarr):
        lc[i] = np.sum(fl[lab == regnum])
        #lc[i] = np.sum(fl[np.where(lab == 1)]
        momim = fl[momlims[0]:momlims[1],
                    momlims[2]:momlims[3]]
        momim[~np.isfinite(momim)] == 0.0
        xbar[i], ybar[i], cov = intertial_axis(momim)


    if return_qual:
        return None
    else:
        return (time,lc, xbar / np.nanmean(xbar), ybar / np.nanmean(xbar), regnum)
       
if __name__=="__main__":

    
    #download K2 TPF from MAST
    path = "."  #Current path
    ar = mastio.K2Archive(path)

    k2id = 206103150 #WASP-47
    campaign = 3

    #extract light curve produced by project office
    fits, fits_hdr = ar.getLongCadence(k2id, campaign, header=True)
    fits_time = fits['TIME']
    fits_flux = fits['SAP_FLUX']
    n_fits_flux = fits_flux/np.median(fits_flux)
    
    #extract flux from TPF
    tpf,tpf_hdr = ar.getLongTpf(k2id, campaign, header=True)
    tpf_arr = getTargetPixelArrayFromFits(tpf, tpf_hdr, column='FLUX')

    
    #measure flux over all finite pixels
    tpf_t_all = tpf['TIME']
    tpf_lc_all = lightCurveFromTpfCube(tpf_arr, mask=None)

    
    #define optimal aperture using a specific set of pixels
    tpf_time = tpf['TIME']
    tpf_flux = tpf['FLUX']
    tpf_quality = tpf['QUALITY']
    time,lc,xbar,ybar,regnum = optimalAperture(tpf_time, tpf_flux, tpf_quality, qual_cut=False, return_qual=False, toss_resat=True, bg_cut=5, skip_cads=None, skip=None)
    tpf_t_opt = time
    tpf_lc_opt = lc

    
    #normalize fluxes
    n_flux_all = tpf_lc_all/np.median(tpf_lc_all)
    n_flux_opt = tpf_lc_opt/np.median(tpf_lc_opt)
    

    goodval_1,= np.where(n_flux_all > 0)
    good_1 = np.isfinite(tpf_lc_opt)
    bad_1 = ~np.isfinite(tpf_lc_opt)
    
    #print "std from optimal pixels"
    #print np.std(n_flux_opt[goodval_1])
    np.savetxt('lc_opt_'+str(k2id)+'_c'+str(campaign)+'.txt',np.c_[tpf_t_opt,tpf_lc_opt,xbar,ybar,bad_1])

    goodval_2,= np.where(n_flux_all > 0)
    good_2 = np.isfinite(tpf_lc_all)
    bad_2 = ~np.isfinite(tpf_lc_all)

    #print "std from all pixels"
    #print np.std(n_flux_all[goodval_2])
    np.savetxt('lc_all_'+str(k2id)+'_c'+str(campaign)+'.txt',np.c_[tpf_t_all,tpf_lc_all,xbar,ybar,bad_2])

    goodval_3,= np.where(n_fits_flux > 0)
    good_3 = np.isfinite(fits_flux)
    bad_3 = ~np.isfinite(fits_flux)

    #print "std from project pixels"
    #print np.std(n_fits_flux[goodval_3])
    np.savetxt('lc_project_'+str(k2id)+'_c'+str(campaign)+'.txt',np.c_[fits_time,fits_flux,xbar,ybar,bad_3])

    
    #plot normalized flux from all pixels versus flux from optimal aperture
    plt.clf()


    #visual comparison of flux from all (red) vs optimal (blue) pixels
    all_flux = np.concatenate([n_fits_flux[goodval_3],n_flux_all[goodval_2],n_flux_opt[goodval_1]])

    plt.plot(fits_time[goodval_3],n_fits_flux[goodval_3],'go',markersize=3.0,label='project pix')
    plt.plot(tpf_t_all[goodval_2],n_flux_all[goodval_2],'ro',markersize=3.0,label='all pix')
    plt.plot(tpf_t_opt[goodval_1],n_flux_opt[goodval_1],'bo',markersize=3.0,label='opt pix')
    plt.xlim(min(tpf_t_all),max(tpf_t_all))
    plt.ylim(min(all_flux),max(all_flux))

    plt.xlabel('BKJD')
    plt.ylabel('Normalized Flux')
    plt.legend(loc='best')

    plt.savefig('lc_full_'+str(k2id)+'_c'+str(campaign)+'.png')
    plt.show()
    plt.clf()
    
    #visual comparison of flux from all (red) vs optimal (blue) pixels
    plt.plot(fits_time,n_fits_flux,'go',markersize=3.0,label='project pix')
    plt.plot(tpf_t_all,n_flux_all,'ro',markersize=3.0,label='all pix')
    plt.plot(tpf_t_opt,n_flux_opt,'bo',markersize=3.0,label='opt pix')
    plt.xlim(min(tpf_t_all),max(tpf_t_all))
    plt.ylim(0.99,1.005)

    plt.xlabel('BKJD')
    plt.ylabel('Normalized Flux')
    plt.legend(loc='best')

    plt.savefig('lc_zoom_'+str(k2id)+'_c'+str(campaign)+'.png')
    plt.show()
    plt.clf()

    figs = plt.figure(1)

    #plot normalized flux from optimal pixels
    plt.subplot(311)
    plt.plot(tpf_t_opt,n_flux_opt,'bo',markersize=3.0,label='opt pix')
    plt.xlim(min(tpf_t_all),max(tpf_t_all))
    plt.ylim(0.96,1.02)
    plt.legend(loc='best')
    plt.ylabel('Normalized Flux')
    
    #plot normalized flux from all pixels
    plt.subplot(312)
    plt.plot(tpf_t_all,n_flux_all,'ro',markersize=3.0,label='all pix')
    plt.xlim(min(tpf_t_all),max(tpf_t_all))
    plt.ylim(0.96,1.02)
    plt.legend(loc='best')
    plt.ylabel('Normalized Flux')

    #plot normalized flux from project office light curve
    plt.subplot(313)
    plt.plot(fits_time,n_fits_flux,'go',markersize=3.0,label='project pix')
    plt.xlim(min(tpf_t_all),max(tpf_t_all))
    plt.ylim(0.96,1.02)
    plt.legend(loc='best')
    plt.ylabel('Normalized Flux')

    plt.xlabel('BKJD')

    plt.savefig('lc_zoom_split_'+str(k2id)+'_c'+str(campaign)+'.png')
    plt.show()
    plt.clf()

    #visual comparison of flux from all (red) vs optimal (blue) pixels
    plt.plot(xbar,ybar,'bo',markersize=3.0)
    plt.xlim(min(xbar),max(xbar))
    plt.ylim(min(ybar),max(ybar))

    plt.xlabel('x centroid')
    plt.ylabel('y centroid')

    plt.savefig('centroid_'+str(k2id)+'_c'+str(campaign)+'.png')
    plt.show()
    plt.clf()
