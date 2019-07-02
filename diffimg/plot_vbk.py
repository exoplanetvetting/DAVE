__version__ = "$Id$"
__URL__ = "$URL$"


import matplotlib.pyplot as mp
import numpy as np
import math

import dave.fileio.kplrfits as kplrfits
import dave.misc.covar as covar
#import diffimg
import dave.diffimg as diffimg
import dave.fileio.mastio as mastio
from dave.diffimg.generate_images import generateImages
from dave.tessPipeline.pertransitcentroids import generateDiffImg
from dave.tessPipeline.pertransitcentroids import getIngressEgressCadences
from dave.tessPipeline.pertransitcentroids import measureCentroidShift
from dave.diffimg.psfFit_vbk import psfCentroids_vbk
import astropy.coordinates as coord
import astropy.units as u
from scipy.signal import medfilt
from numpy.matlib import repmat
from matplotlib import gridspec
from astropy import coordinates, units as u, wcs
from astroquery.skyview import SkyView
from astroquery.vizier import Vizier
import astropy.units as u
from scipy.ndimage import rotate
#from reproject import reproject_interp, reproject_exact, reproject_to_healpix, reproject_from_healpix
from astropy.wcs import WCS
#import pywcsgrid2


from matplotlib.colors import LinearSegmentedColormap

cm_data = [[0.2081, 0.1663, 0.5292], [0.2116238095, 0.1897809524, 0.5776761905], 
 [0.212252381, 0.2137714286, 0.6269714286], [0.2081, 0.2386, 0.6770857143], 
 [0.1959047619, 0.2644571429, 0.7279], [0.1707285714, 0.2919380952, 
  0.779247619], [0.1252714286, 0.3242428571, 0.8302714286], 
 [0.0591333333, 0.3598333333, 0.8683333333], [0.0116952381, 0.3875095238, 
  0.8819571429], [0.0059571429, 0.4086142857, 0.8828428571], 
 [0.0165142857, 0.4266, 0.8786333333], [0.032852381, 0.4430428571, 
  0.8719571429], [0.0498142857, 0.4585714286, 0.8640571429], 
 [0.0629333333, 0.4736904762, 0.8554380952], [0.0722666667, 0.4886666667, 
  0.8467], [0.0779428571, 0.5039857143, 0.8383714286], 
 [0.079347619, 0.5200238095, 0.8311809524], [0.0749428571, 0.5375428571, 
  0.8262714286], [0.0640571429, 0.5569857143, 0.8239571429], 
 [0.0487714286, 0.5772238095, 0.8228285714], [0.0343428571, 0.5965809524, 
  0.819852381], [0.0265, 0.6137, 0.8135], [0.0238904762, 0.6286619048, 
  0.8037619048], [0.0230904762, 0.6417857143, 0.7912666667], 
 [0.0227714286, 0.6534857143, 0.7767571429], [0.0266619048, 0.6641952381, 
  0.7607190476], [0.0383714286, 0.6742714286, 0.743552381], 
 [0.0589714286, 0.6837571429, 0.7253857143], 
 [0.0843, 0.6928333333, 0.7061666667], [0.1132952381, 0.7015, 0.6858571429], 
 [0.1452714286, 0.7097571429, 0.6646285714], [0.1801333333, 0.7176571429, 
  0.6424333333], [0.2178285714, 0.7250428571, 0.6192619048], 
 [0.2586428571, 0.7317142857, 0.5954285714], [0.3021714286, 0.7376047619, 
  0.5711857143], [0.3481666667, 0.7424333333, 0.5472666667], 
 [0.3952571429, 0.7459, 0.5244428571], [0.4420095238, 0.7480809524, 
  0.5033142857], [0.4871238095, 0.7490619048, 0.4839761905], 
 [0.5300285714, 0.7491142857, 0.4661142857], [0.5708571429, 0.7485190476, 
  0.4493904762], [0.609852381, 0.7473142857, 0.4336857143], 
 [0.6473, 0.7456, 0.4188], [0.6834190476, 0.7434761905, 0.4044333333], 
 [0.7184095238, 0.7411333333, 0.3904761905], 
 [0.7524857143, 0.7384, 0.3768142857], [0.7858428571, 0.7355666667, 
  0.3632714286], [0.8185047619, 0.7327333333, 0.3497904762], 
 [0.8506571429, 0.7299, 0.3360285714], [0.8824333333, 0.7274333333, 0.3217], 
 [0.9139333333, 0.7257857143, 0.3062761905], [0.9449571429, 0.7261142857, 
  0.2886428571], [0.9738952381, 0.7313952381, 0.266647619], 
 [0.9937714286, 0.7454571429, 0.240347619], [0.9990428571, 0.7653142857, 
  0.2164142857], [0.9955333333, 0.7860571429, 0.196652381], 
 [0.988, 0.8066, 0.1793666667], [0.9788571429, 0.8271428571, 0.1633142857], 
 [0.9697, 0.8481380952, 0.147452381], [0.9625857143, 0.8705142857, 0.1309], 
 [0.9588714286, 0.8949, 0.1132428571], [0.9598238095, 0.9218333333, 
  0.0948380952], [0.9661, 0.9514428571, 0.0755333333], 
 [0.9763, 0.9831, 0.0538]]

parula_map = LinearSegmentedColormap.from_list('parula', cm_data)

def plotWrapper(clip):
    """Wrapper function for difference image centroid diagnostic plots

    Call this function from the exporter.

    Inputs:
    -----------
    clip
        A clipboard object. Should have the following keys: serve, detrend, diffImg, rollPhase, trapFit

    Returns:
    -----------
    Two figure handles. The zeroth figure handle shows the flux and rollPhase
    plot, the first shows the centroid offset plot

    Outputs:
    ------------
    Two figures are produced.

    """
    time = clip['serve.time']
    qFlags = clip['serve.flags']

    flux = clip['detrend.flux_frac']
    flags = clip['detrend.flags']

    period_days = clip['trapFit.period_days']
    epoch_bkjd = clip['trapFit.epoch_bkjd']
    duration_hrs = clip['trapFit.duration_hrs']
    epic = clip['value']

    cube = clip['serve.cube']

    inTransitIndices = kplrfits.markTransitCadences(time, period_days, epoch_bkjd, duration_hrs/24., flags=flags)

#    tce = clip['eventList'][0]
#    period = tce['trapFit.period_days']
#    epoch = tce['trapFit.epoch_bkjd']
#    duration_hrs = tce['trapFit.duration_hrs']

    if clip['config.detrendType'] != "tess_2min" and (clip['config.detrendType'] != "eleanor"):
        rollPhase = clip['rollPhase.rollPhase']
        centroids = clip['diffImg.centroid_timeseries']
        goodCentroidIndices = centroids[ centroids[:,1]>1, 0].asarray().astype(int)

        fig1 = mp.figure(1)
        mp.clf()
        multiPanelPlotDiffImgCentroidsDiagnostic(time, flux, flags, rollPhase, inTransitIndices, goodCentroidIndices, qFlags)

        fig2 = mp.figure(2)
        mp.clf()
        try:
            titleStr = PLOT_CENTROID_OFFSETS_VBK(clip)
#        titleStr = "EPIC: %i  %s" %(epic, titleStr)
        except ValueError as e:
            titleStr = "Error: %s" %(e)

    elif clip['config.detrendType'] == "tess_2min":
        idx_arr_tmp_ = np.asarray(inTransitIndices)
        goodCentroidIndices = np.linspace(0,len(idx_arr_tmp_), len(idx_arr_tmp_)+1)
        goodCentroidIndices = goodCentroidIndices[:-1]

        fig1 = mp.figure(1)
        mp.clf()
        xx = PLOT_CENTROIDS_TESS(clip)#PLOT_DIFF_IMG_TESS(clip)#

        fig2 = mp.figure(2)
        mp.clf()
        xx = PLOT_INDIV_IMG_TESS(clip)

    return fig1, fig2
##
##
##
def PLOT_CENTROID_OFFSETS_VBK(clip):

    from astropy import coordinates, units as u, wcs
    from astroquery.skyview import SkyView
    from astroquery.vizier import Vizier
    import astropy.units as u
    import math
    from scipy.ndimage import rotate
#    from reproject import reproject_interp, reproject_exact, reproject_to_healpix, reproject_from_healpix
    from astropy.wcs import WCS
#    import pywcsgrid2

    time = clip['serve.time']
    qFlags = clip['serve.flags']

    flux = clip['detrend.flux_frac']
    flags = clip['detrend.flags']

    centroids = clip['diffImg.centroid_timeseries']
#    rollPhase = clip['rollPhase.rollPhase']
    period_days = clip['trapFit.period_days']
    epoch_bkjd = clip['trapFit.epoch_bkjd']
    duration_hrs = clip['trapFit.duration_hrs']
    epic = clip['value']

    cube = clip['serve.cube']

    hdr_ = clip['serve.tpfHeader']
    col_zero_, row_zero_ = int(hdr_['1CRV4P']), int(hdr_['2CRV4P'])

    epic_Col, epic_Row = col_zero_ + int(hdr_['1CRPX4']), row_zero_ + int(hdr_['2CRPX4'])

    def k2_ConvertHeaderWCS(tpf_header):
        funny_keywords = {'1CTYP4': 'CTYPE1',
                        '2CTYP4': 'CTYPE2',
                        '1CRPX4': 'CRPIX1',
                        '2CRPX4': 'CRPIX2',
                        '1CRVL4': 'CRVAL1',
                        '2CRVL4': 'CRVAL2',
                        '1CUNI4': 'CUNIT1',
                        '2CUNI4': 'CUNIT2',
                        '1CDLT4': 'CDELT1',
                        '2CDLT4': 'CDELT2',
                        '11PC4': 'PC1_1',
                        '12PC4': 'PC1_2',
                        '21PC4': 'PC2_1',
                        '22PC4': 'PC2_2'}
        mywcs = {}
        for oldkey, newkey in funny_keywords.items():
                mywcs[newkey] = tpf_header[oldkey]

        return wcs.WCS(mywcs)

    mywcs_ = k2_ConvertHeaderWCS(hdr_)
#
#
    inTransitIndices = kplrfits.markTransitCadences(time, period_days, epoch_bkjd, duration_hrs/24., flags=flags)

    oot_cadence_ = np.where((inTransitIndices == False) & (flags == False))
    oot_mean_img_ = np.nanmean(cube[oot_cadence_], axis=0)

    itr_cadence_ = np.where(inTransitIndices == True)
    itr_mean_img_ = np.nanmean(cube[itr_cadence_], axis=0)

    diff_mean_img_ = oot_mean_img_ - itr_mean_img_

    ss_ = oot_mean_img_.shape

#    disp = lambda x: mp.imshow(x, cmap=mp.cm.binary, origin = "bottom", interpolation="nearest")
    extent_ = [col_zero_, col_zero_ + ss_[1], row_zero_, row_zero_ + ss_[0]]
#    cmap_ = lambda x: mp.get_cmap('binary', int(np.max(data))-int(np.min(data))+1)
    disp = lambda x: mp.imshow(x, cmap=mp.get_cmap('binary', 512), origin = "bottom", interpolation="nearest", extent = extent_)
#
# GET CENTROIDS

    idx = centroids[:,1] > 0
    cin = centroids[idx, 0]

    ootCol = centroids[idx, 1]# - col_zero_# - 1
    ootRow = centroids[idx, 2]# - row_zero_# - 1

    #itr => in transit
    diffCol = centroids[idx, 3]# - col_zero_# - 1 
    diffRow = centroids[idx, 4]# - row_zero_# - 1

    diffC = (ootCol - diffCol)# + np.median(diffCol)
    diffR = (ootRow - diffRow)# + np.median(diffRow)

    itrCol, itrRow = diffCol, diffRow

#
    xmin_ = np.min(np.hstack((ootCol, diffCol)))
    xmax_ = np.max(np.hstack((ootCol, diffCol)))

    ymin_ = np.min(np.hstack((ootRow, diffRow)))
    ymax_ = np.max(np.hstack((ootRow, diffRow)))

# START PLOTTING
#
    ax1 = mp.subplot(221)

    disp(oot_mean_img_)

#    mp.scatter(diffC, diffR, marker='o', c=cin, s=64, linewidths=0, cmap=mp.cm.RdYlBu)
#    mp.plot(diffC, diffR, 'ro', ms= 6)

#    mp.axhline(0, color='k', lw=.6)
#    mp.axvline(0, color='k', lw=.5)
#    mp.plot(0.,0., '*', ms=40, color='yellow')
    ax1.plot(ootCol, ootRow, 'c*', ms=8)#, mec = 'm')#, color='yellow')
    ax1.plot(np.mean(ootCol), np.mean(ootRow), 'c*', ms=14, label = 'AVG_OOT')#, mec = 'm')

    ax1.plot(itrCol, itrRow, 'mo', ms=3)#, mec = 'c')
    ax1.plot(np.mean(itrCol), np.mean(itrRow), 'mo', ms=6, label = 'AVG_DIFF')#, mec = 'c')
#
    covar.plotErrorEllipse(ootCol, ootRow, color='c', ms=14, marker = '*', mfc = 'c')#, mec = 'm')
    covar.plotErrorEllipse(itrCol, itrRow, color='m', ms=14, marker = 'o', mfc = 'm')#, mec = 'c')
#    covar.plotErrorEllipse(diffC, diffR, color='#888888', ms=3)

    ax1.plot(epic_Col, epic_Row, 'xy', mew=3, ms = 10, label = 'EPIC')

    mp.xlabel(r"$\Delta$ Column (pixels)")
    mp.ylabel(r"$\Delta$ Row (pixels)")

    mp.legend(loc = 'best', fontsize = 8)

#    xmin_ = np.min(np.asarray([ootCol, diffCol]))
#    xmax_ = np.max(np.asarray([ootCol, diffCol]))

#    ymin_ = np.min(np.asarray([ootRow, diffRow]))
#    ymax_ = np.max(np.asarray([ootRow, diffRow]))

#    mp.xlim(0.8*xmin_, 1.2*xmax_)
#    mp.ylim(0.8*ymin_, 1.2*ymax_)

#    probOffset, chiSq = covar.computeProbabilityOfObservedOffset(diffC, diffR)
#    titleStr = "Prob. On Target: %.1e: $\chi^2$: %.3f" %(1-probOffset, chiSq)
#    cb = mp.colorbar()
#    cb.set_label("Time (BKJD)")

    #Ensure some padding around the origin so the symbol is never
    #at edge of the plot
#    axl = list(mp.axis())
#    axl[0] = min(axl[0], -0.1)
#    axl[1] = max(axl[1], +0.1)
#    axl[2] = min(axl[2], -0.1)
#    axl[3] = max(axl[3], +0.1)
#    mp.axis(axl)

    mp.subplot(222)
##    titleStr = plotCentroidOffsets(centroids)
##
##
#    mp.scatter(diffC, diffR, marker='o', c=cin, s=64, linewidths=0, cmap=mp.cm.RdYlBu)
#    mp.plot(ootCol, ootRow, 'c*', ms = 10, mec = 'k')#cmap=mp.cm.Wistia)#spring)
#    mp.plot(itrCol, itrRow, 'mo', ms = 10, mec = 'k')#cmap=mp.cm.winter)

    mp.plot(np.mean(ootCol),np.mean(ootRow), 'c*', ms=20, label = 'OOT')
    mp.plot(np.mean(itrCol),np.mean(itrRow), 'mo', ms=20, label = 'DIFF')
    mp.scatter(ootCol, ootRow, marker='*', c=cin, s=64, linewidths=0, cmap=mp.cm.RdYlBu)
    mp.scatter(itrCol, itrRow, marker='o', c=cin, s=64, linewidths=0, cmap=mp.cm.RdYlBu)

    cb = mp.colorbar()
    cb.set_label("Cadence")

#    mp.axhline(0, color='k', lw=.5)
#    mp.axvline(0, color='k', lw=.5)

#    covar.plotErrorEllipse(diffC, diffR, color='#888888', ms=20)
    covar.plotErrorEllipse(ootCol, ootRow, color='c', ms=20, marker = '*', mfc = 'c')
    covar.plotErrorEllipse(itrCol, itrRow, color='m', ms=20, marker = 'o', mfc = 'm')

    mp.xlabel(r"$\Delta$ Column (pixels)")

#    mp.ylabel(r"$\Delta$ Row (pixels)")

    probOffset, chiSq = covar.computeProbabilityOfObservedOffset(diffC, diffR)
    titleStr = "Prob. On Target: %.1e: $\chi^2$: %.3f" %(1-probOffset, chiSq)

    mp.legend(loc = 'best')
#    cb = mp.colorbar()
#    cb.set_label("Time (BKJD)")

    mp.xlim(xmin_-0.2, xmax_+0.2)
    mp.ylim(ymin_-0.2, ymax_+0.2)
    mp.tight_layout()
#    mp.locator_params(axis='y', nbins=4)
#    mp.locator_params(axis='x', nbins=4)


    #Ensure some padding around the origin so the symbol is never
    #at edge of the plot
#    axl = list(mp.axis())
#    axl[0] = min(axl[0], -0.1)
#    axl[1] = max(axl[1], +0.1)
#    axl[2] = min(axl[2], -0.1)
#    axl[3] = max(axl[3], +0.1)
#    mp.axis(axl)
#
#

    try:
        ax3 = mp.subplot(223, projection = mywcs_)#mywcs_)

        ra_, dec_ = hdr_['RA_OBJ'], hdr_['DEC_OBJ']
        center_ = coordinates.SkyCoord(ra_, dec_, unit=(u.deg, u.deg), frame='icrs')

        img_survey = SkyView.get_images(position=center_, survey='2MASS-J', radius=1*u.arcmin)
        pix_survey = img_survey[0][0].data
        hdr_survey = img_survey[0][0].header

#    levels_ = np.linspace(np.min(pix_survey),np.percentile(pix_survey,95),40)
#    levels_ = [0.9, 0.5, 0.1, 0.0]

        inverted_pix_survey = np.max(pix_survey) - pix_survey
        inverted_pix_survey = pix_survey#inverted_pix_survey/np.max(inverted_pix_survey)

        levels_ = np.linspace(np.min(inverted_pix_survey),np.percentile(inverted_pix_survey,99),10)
#    levels_ = np.logspace(np.min(inverted_pix_survey),np.max(inverted_pix_survey),10)
#    levels_ = [0.0,0.0001,0.0005,0.001,0.005,0.01,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95,0.99,1.]

#    array, footprint = reproject_interp(pix_survey, mywcs_, shape_out = pix_survey.shape)#catalog_img_.shape)#ss_)

#    disp(pix_survey)#.T)
#    ax3.imshow(pix_survey, cmap=mp.get_cmap('binary',512))#, transform = wcs.WCS(hdr_survey))
        ax3.contourf(inverted_pix_survey, transform=ax3.get_transform(wcs.WCS(hdr_survey)), levels = levels_, cmap=mp.get_cmap('binary',256))
#    ax.scatter(ra_, dec_, 'xy', mew=3, ms = 10, label = 'EPIC', transform=ax.get_transform(wcs.WCS(hdr_survey)))
#    print ax3.get_xlim(), ax3.get_ylim()
#
#
        mp.tight_layout()
    except:
        mp.subplot(223)

    mp.subplot(224)
    titleStr_ = plotCentroidOffsets(centroids)
#    disp(diff_mean_img_)
#    mp.scatter(diffC, diffR, marker='s', c=cin, s=64, linewidths=0, cmap=mp.cm.RdYlBu)
#    covar.plotErrorEllipse(diffC, diffR, color='#888888', ms=20, marker = 's')

#    mp.scatter(ootCol, ootRow, marker='*', c=cin, s=64, linewidths=0, cmap=mp.cm.Wistia)#spring)
#    covar.plotErrorEllipse(ootCol, ootRow, color='#888888', ms=20, marker = '*', mfc = 'm')

#    mp.scatter(itrCol, itrRow, marker='o', c=cin, s=64, linewidths=0, cmap=mp.cm.winter) 
#    covar.plotErrorEllipse(itrCol, itrRow, color='#888888', ms=20, marker = 'o', mfc = 'b')

    return titleStr
##
##
##
def PLOT_CENTROIDS_TESS(clip):

    out = clip['diffImgCentroids']
    centroids = clip['diffImgCentroids.results']

    ootCol_prf, ootRow_prf = np.mean([centroids[:,0],centroids[:,4]], axis = 0), np.mean([centroids[:,1],centroids[:,5]], axis = 0)
    diffCol_prf, diffRow_prf = centroids[:,2], centroids[:,3]

#    print(ootCol_prf, ootRow_prf)
#    print(diffCol_prf, diffRow_prf)

#    itrCol_psf, itrRow_psf, ootCol_psf, ootRow_psf, diffCol_psf, diffRow_psf = psfCentroids_vbk(clip)
#    diffCol_prf, diffRow_prf = diffCol_psf, diffRow_psf
#    ootCol_prf, ootRow_prf = ootCol_psf, ootRow_psf

    diffC = (ootCol_prf - diffCol_prf)
    diffR = (ootRow_prf - diffRow_prf)

#    itr_mean_img_, oot_mean_img_, diff_mean_img_, itr_mean_cube_, oot_mean_cube_, diff_mean_cube_, transit_number_ = generateImages(clip)
#    itr_mean_img_, oot_mean_img_, diff_mean_img_ = np.asarray(itr_mean_img_), np.asarray(oot_mean_img_), np.asarray(diff_mean_img_)

    cube = clip['serve.cube']
    time = clip['serve.time']
    period_days = clip['trapFit.period_days']
    epoch_days = clip['trapFit.epoch_bkjd']
    duration_days = clip['trapFit.duration_hrs']/24.

    isnan = np.isnan(time)
    time = time[~isnan]
    cube = cube[~isnan]

    transits = getIngressEgressCadences(time, period_days, epoch_days, duration_days)

    oot_img_ = np.zeros((len(transits), cube.shape[1], cube.shape[2]))
    diff_img_ = np.zeros((len(transits), cube.shape[1], cube.shape[2]))

    transit_flags_ = np.zeros(len(transits), dtype = bool)
    for jj in range(len(transits)):
        key = 'transit-%04i' %(jj)
        if out[key]['errorCode'] < 7:
            cin = transits[jj]
            plot = False
            before, after, diff = generateDiffImg(cube, cin, plot=plot)
            diff_img_[jj,:,:] = diff
            oot_img_[jj,:,:] = 0.5*(before + after)
        else:
            transit_flags_[jj] = True

    transits = transits[~transit_flags_]
    diff_img_, oot_img_ = diff_img_[~transit_flags_], oot_img_[~transit_flags_]
    ootCol_prf, ootRow_prf = ootCol_prf[~transit_flags_], ootRow_prf[~transit_flags_]
    diffCol_prf, diffRow_prf = diffCol_prf[~transit_flags_], diffRow_prf[~transit_flags_]
    diffC, diffR = diffC[~transit_flags_], diffR[~transit_flags_]

#    print(oot_img_.shape, oot_img_.shape)
    oot_mean_img_ = np.nanmean(oot_img_, axis = 0)
    diff_mean_img_ = np.nanmean(diff_img_, axis = 0)
    itr_mean_img_ = oot_mean_img_

#    print(diff_img_.shape, oot_img_.shape, oot_mean_img_.shape, diff_mean_img_.shape)
#    xxxx

    hdr_ = clip['serve.tpfHeader']    

    if clip['config.detrendType'] == "eleanor":
    	col_zero_, row_zero_ = int(hdr_['CRPIX1']), int(hdr_['CRPIX2'])
    	epic_Col, epic_Row = col_zero_ + int(hdr_['TPF_W']), row_zero_ + int(hdr_['TPF_H'])

    ss_ = itr_mean_img_.shape

    extent_ = [0, ss_[1],0, ss_[0]]
    disp = lambda x: mp.imshow(x, cmap = parula_map, origin = "bottom", interpolation = "nearest", extent = extent_)

    mp.subplot(221)

    disp(diff_mean_img_)

#    mp.plot(np.mean(diffCol_prf), np.mean(diffRow_prf), 'r*', ms = 12, mec = 'k', label="DIFF")
    mp.plot(diffCol_prf, diffRow_prf, 'r*')#, label="Indiv Diff")
#    mp.plot(np.mean(ootCol_prf), np.mean(ootRow_prf), 'mo', label="OOT")
   
#    for i in range(len(centroids)):
#        mp.text(centroids[i,2], centroids[i,3], ' %i' %(i))

#    covar.plotErrorEllipse(diffCol_prf, diffRow_prf, color='r', ms=14, marker = '*', mfc = 'r')
    probOffset, chiSq = covar.computeProbabilityOfObservedOffset(diffC, diffR)
    titleStr = "Prob. On Target: %.1e; $\chi^2$: %.3f" %(1-probOffset, chiSq)

#    mp.xlabel("Column")
#    mp.ylabel("Row")
#    mp.legend(loc = 'best', fontsize = 10)
#    mp.title('Diff. Image; ' + titleStr, fontsize = 10)
    mp.title('Difference Image')
    mp.colorbar()

    mp.subplot(222)

    disp(oot_mean_img_)
    mp.plot(ootCol_prf, ootRow_prf, 'mo', label="OOT")

#    mp.xlabel("Column")
#    mp.ylabel("Row")
#    mp.legend(loc = 'best', fontsize = 10)
#    mp.title('OOT', fontsize = 10)
    mp.title('Out-of-transit Image')
    mp.colorbar()

    mp.subplot(223)

    disp(itr_mean_img_)
    mp.plot(ootCol_prf, ootRow_prf, 'mo', label="OOT")
#    mp.plot(centroids[:,0], centroids[:,1], 'ko', label="Before")
#    mp.plot(centroids[:,2], centroids[:,3], 'ro', label="Difference")
#    mp.plot(centroids[:,4], centroids[:,5], 'wo', label="After")
#    mp.plot(np.mean(centroids[:,0]), np.mean(centroids[:,1]), 'ks', label="Before")
#    mp.plot(np.mean(centroids[:,4]), np.mean(centroids[:,5]), 'wo', label="After")
#    mp.plot(diffCol+0*col_zero_,diffRow+0*row_zero_,'r*')#,ms=6)
  
#    mp.xlabel("Column")
#    mp.ylabel("Row")
#    mp.legend(loc = 'best')
    mp.title('In-transit Image')
    mp.colorbar()
#
#
#
#

    def k2_ConvertHeaderWCS(tpf_header):
        funny_keywords = {'1CTYP4': 'CTYPE1',
            '2CTYP4': 'CTYPE2', '1CRPX4': 'CRPIX1', '2CRPX4': 'CRPIX2','1CRVL4': 'CRVAL1','2CRVL4': 'CRVAL2','1CUNI4': 'CUNIT1','2CUNI4': 'CUNIT2',
            '1CDLT4': 'CDELT1','2CDLT4': 'CDELT2','11PC4': 'PC1_1','12PC4': 'PC1_2','21PC4': 'PC2_1','22PC4': 'PC2_2'}
        mywcs = {}
        for oldkey, newkey in funny_keywords.items():
            mywcs[newkey] = tpf_header[oldkey]
        return wcs.WCS(mywcs)

    mywcs_ = k2_ConvertHeaderWCS(hdr_)

    
    plot_ = 'SNR'#'img'
    if plot_ == 'img':

    	ax3 = mp.subplot(224, projection = mywcs_)#mywcs_)

    	hdr_ = clip['serve.tpfHeader']
    	ra_, dec_ = hdr_['RA_OBJ'], hdr_['DEC_OBJ']
    	center_ = coordinates.SkyCoord(ra_, dec_, unit=(u.deg, u.deg), frame='icrs')

    	img_survey = SkyView.get_images(position=center_, survey='DSS', radius=3.5*u.arcmin)
    	# NEED TO DO pip install urllib3==1.22 TO GET SOMETHING CORRECT!!!!
    	pix_survey = img_survey[0][0].data
    	hdr_survey = img_survey[0][0].header

    	inverted_pix_survey = np.max(pix_survey) - pix_survey
    	inverted_pix_survey = pix_survey#inverted_pix_survey/np.max(inverted_pix_survey)

    	levels_ = np.linspace(np.min(inverted_pix_survey),np.percentile(inverted_pix_survey,95),10)

    	ax3.contourf(inverted_pix_survey,transform=ax3.get_transform(wcs.WCS(hdr_survey)),levels = levels_,cmap=mp.get_cmap('binary',256))

    elif plot_ == 'SNR':

    	mp.subplot(224)
    	tmp_diff_img_ = 0.01 + diff_mean_img_ - np.min(diff_mean_img_.flatten())
    	idx = np.isfinite(tmp_diff_img_)
    	snr = np.empty_like(tmp_diff_img_)
    	snr[idx] = tmp_diff_img_[idx]/np.sqrt(tmp_diff_img_[idx])
    	disp(snr)
    	mp.colorbar()
    	mp.title("Difference SNR")
#    mp.xlabel("Column")


    return titleStr
##
##
##
def PLOT_DIFF_IMG_TESS(clip):

    out = clip['diffImgCentroids']
    centroids = clip['diffImgCentroids.results']
    ootCol_prf, ootRow_prf = np.mean([centroids[:,0],centroids[:,4]], axis = 0), np.mean([centroids[:,1],centroids[:,5]], axis = 0)
    diffCol_prf, diffRow_prf = centroids[:,2], centroids[:,3]

    itrCol, itrRow = ootCol_prf - diffCol_prf, ootRow_prf - diffRow_prf
    ootCol, ootRow, diffCol, diffRow = ootCol_prf, ootRow_prf, diffCol_prf, diffRow_prf

#    itrCol, itrRow, ootCol, ootRow, diffCol, diffRow = psfCentroids_vbk(clip)

    diffC = (ootCol - diffCol)
    diffR = (ootRow - diffRow)

#    itr_mean_img_, oot_mean_img_, diff_mean_img_, itr_mean_cube_, oot_mean_cube_, diff_mean_cube_, transit_number_ = generateImages(clip)
    itr_mean_img_, oot_mean_img_, diff_mean_img_ = generateImages(clip)
    itr_mean_img_, oot_mean_img_, diff_mean_img_ = np.asarray(itr_mean_img_), np.asarray(oot_mean_img_), np.asarray(diff_mean_img_)

    hdr_ = clip['serve.tpfHeader']    

    if clip['config.detrendType'] == "tess_2min":
    	col_zero_, row_zero_ = int(hdr_['1CRV4P']), int(hdr_['2CRV4P'])
    	epic_Col, epic_Row = col_zero_ + int(hdr_['1CRPX4']), row_zero_ + int(hdr_['2CRPX4'])

    if clip['config.detrendType'] == "eleanor":
    	col_zero_, row_zero_ = int(hdr_['CRPIX1']), int(hdr_['CRPIX2'])
    	epic_Col, epic_Row = col_zero_ + int(hdr_['TPF_W']), row_zero_ + int(hdr_['TPF_H'])
#
#
# START PLOTTING CENTROIDS
#
#
    ss_ = oot_mean_img_.shape
#    extent_ = [col_zero_, col_zero_ + ss_[1], row_zero_, row_zero_ + ss_[0]]
#    disp = lambda x: mp.imshow(x, cmap=mp.cm.YlGnBu_r, origin = "bottom", interpolation = "nearest", extent = extent_)
    extent_ = [0, ss_[1],0, ss_[0]]
    disp = lambda x: mp.imshow(x, cmap= parula_map, origin = "bottom", interpolation = "nearest", extent = extent_)

    mp.subplot(221)

    plot_ = 'plot_images_'#'plot_flux_'#

    if plot_ == 'plot_flux_':
        mp.plot(time[oot_cadence_], flux[oot_cadence_], 'k.')
        mp.plot(time[itr_cadence_], flux[itr_cadence_], 'r.')
        mp.xlim(1356., 1357.)
    elif plot_ == 'plot_images_':
        disp(diff_mean_img_)
        mp.plot(diffCol,diffRow,'r*',ms=6)
        mp.colorbar()
        mp.clim(vmin=0)
        covar.plotErrorEllipse(diffCol, diffRow, color='r', ms=14, marker = '*', mfc = 'r')#, mec = 'c')
        probOffset, chiSq = covar.computeProbabilityOfObservedOffset(diffC, diffR)
        titleStr = "Prob. On Target: %.1e; $\chi^2$: %.3f" %(1-probOffset, chiSq)

    mp.title("Difference Image")#: cent @ ["+str(round(mean_diffRow +1*row_zero_,2))+','+str(round(mean_diffCol +1*col_zero_,2))+']')

    mp.subplot(222)
    disp(oot_mean_img_)
    mp.plot(ootCol,ootRow,'mo',ms=6)# if use_before_after_images_per_transit_=='n' else mp.plot(mean_ootCol,mean_ootRow,'mo',ms=9)
#    mp.plot(itrCol, itrRow, 'mo', ms = 4)
#    mp.plot(diffCol, diffRow, 'gs', ms = 4)
    mp.colorbar()
    mp.title("OOT Image")#: cent @ ["+str(round(mean_ootRow +1*row_zero_,2))+','+str(round(mean_ootCol + 1*col_zero_,2)) + ']')
#
#
#
#    
    plt_what_ = 'itr'#'archive_img'#

    if plt_what_ != 'archive_img':
        mp.subplot(223)

        disp(itr_mean_img_)
#    mp.plot(ootCol, ootRow, 'c*', ms = 4)
        mp.plot(itrCol,itrRow,'mo',ms=6)# if use_before_after_images_per_transit_=='n' else mp.plot(mean_itrCol,mean_itrRow,'gx',ms=9)
#    mp.plot(diffCol, diffRow, 'gs', ms = 4)
        mp.colorbar()
        mp.title("In-transit Image")#: cent @ ["+str(round(mean_itrRow + 1*row_zero_,2))+','+str(round(mean_itrCol + 1*col_zero_,2)) + ']')

    else:

    	try:
            from astropy import coordinates, units as u, wcs
            from astroquery.skyview import SkyView
            from astroquery.vizier import Vizier
            import astropy.units as u
            import math
            from scipy.ndimage import rotate
#            from reproject import reproject_interp, reproject_exact, reproject_to_healpix, reproject_from_healpix
            from astropy.wcs import WCS
#    		import pywcsgrid2

            hdr_ = clip['serve.tpfHeader']
            col_zero_, row_zero_ = int(hdr_['1CRV4P']), int(hdr_['2CRV4P'])
            epic_Col, epic_Row = col_zero_ + int(hdr_['1CRPX4']), row_zero_ + int(hdr_['2CRPX4'])

            def k2_ConvertHeaderWCS(tpf_header):
                funny_keywords = {'1CTYP4': 'CTYPE1',
                    '2CTYP4': 'CTYPE2',
		    '1CRPX4': 'CRPIX1',
                    '2CRPX4': 'CRPIX2',
                    '1CRVL4': 'CRVAL1',
                    '2CRVL4': 'CRVAL2',
                    '1CUNI4': 'CUNIT1',
                    '2CUNI4': 'CUNIT2',
                    '1CDLT4': 'CDELT1',
                    '2CDLT4': 'CDELT2',
                    '11PC4': 'PC1_1',
                    '12PC4': 'PC1_2',
                    '21PC4': 'PC2_1',
                    '22PC4': 'PC2_2'}
                mywcs = {}
                for oldkey, newkey in funny_keywords.items():
                    mywcs[newkey] = tpf_header[oldkey]

                return wcs.WCS(mywcs)

            mywcs_ = k2_ConvertHeaderWCS(hdr_)

            ax3 = mp.subplot(223, projection = mywcs_)#mywcs_)

            ra_, dec_ = hdr_['RA_OBJ'], hdr_['DEC_OBJ']

            print(ra_, dec_)
            center_ = coordinates.SkyCoord(ra_, dec_, unit=(u.deg, u.deg), frame='icrs')

            img_survey = SkyView.get_images(position=center_, survey='2MASS-J', radius=3.5*u.arcmin)
            pix_survey = img_survey[0][0].data
            hdr_survey = img_survey[0][0].header

            inverted_pix_survey = np.max(pix_survey) - pix_survey
            inverted_pix_survey = pix_survey#inverted_pix_survey/np.max(inverted_pix_survey)

            levels_ = np.linspace(np.min(inverted_pix_survey),np.percentile(inverted_pix_survey,99),10)

            ax3.contourf(inverted_pix_survey,transform=ax3.get_transform(wcs.WCS(hdr_survey)),levels = levels_,cmap=mp.get_cmap('binary',256))
    	except:
        	mp.subplot(223)
#
#
#
#
    plt_what_ = 'archive_img'#'snr'#
    if plt_what_ != 'archive_img':
        mp.subplot(224)

        tmp_diff_img_ = 0.01 + diff_mean_img_ - np.min(diff_mean_img_.flatten())
        idx = np.isfinite(tmp_diff_img_)
        snr = np.empty_like(tmp_diff_img_)
        snr[idx] = tmp_diff_img_[idx]/np.sqrt(tmp_diff_img_[idx])
        disp(snr)
        mp.colorbar()
        mp.title("Difference SNR")
    elif plt_what_ == 'archive_img':
        try:
            from astropy import coordinates, units as u, wcs
            from astroquery.skyview import SkyView
            from astroquery.vizier import Vizier
#    		import astropy.units as u
            import math
            from scipy.ndimage import rotate
            from reproject import reproject_interp, reproject_exact, reproject_to_healpix, reproject_from_healpix
            from astropy.wcs import WCS
#    		import pywcsgrid2

            hdr_ = clip['serve.tpfHeader']
            col_zero_, row_zero_ = int(hdr_['1CRV4P']), int(hdr_['2CRV4P'])
            epic_Col, epic_Row = col_zero_ + int(hdr_['1CRPX4']), row_zero_ + int(hdr_['2CRPX4'])

            def k2_ConvertHeaderWCS(tpf_header):
                funny_keywords = {'1CTYP4': 'CTYPE1',
                        '2CTYP4': 'CTYPE2',
                        '1CRPX4': 'CRPIX1',
                        '2CRPX4': 'CRPIX2',
                        '1CRVL4': 'CRVAL1',
                        '2CRVL4': 'CRVAL2',
                        '1CUNI4': 'CUNIT1',
                        '2CUNI4': 'CUNIT2',
                        '1CDLT4': 'CDELT1',
                        '2CDLT4': 'CDELT2',
                        '11PC4': 'PC1_1',
                        '12PC4': 'PC1_2',
                        '21PC4': 'PC2_1',
                        '22PC4': 'PC2_2'}
                mywcs = {}
                for oldkey, newkey in funny_keywords.items():
                    mywcs[newkey] = tpf_header[oldkey]

                return wcs.WCS(mywcs)

            mywcs_ = k2_ConvertHeaderWCS(hdr_)

            ax3 = mp.subplot(224, projection = mywcs_)#mywcs_)

            ra_, dec_ = hdr_['RA_OBJ'], hdr_['DEC_OBJ']

#		print ra_, dec_
            center_ = coordinates.SkyCoord(ra_, dec_, unit=(u.deg, u.deg), frame='icrs')

            img_survey = SkyView.get_images(position=center_, survey='2MASS-J', radius=3.5*u.arcmin)
		# NEED TO DO pip install urllib3==1.22 TO GET SOMETHING CORRECT!!!!
            pix_survey = img_survey[0][0].data
            hdr_survey = img_survey[0][0].header

            inverted_pix_survey = np.max(pix_survey) - pix_survey
            inverted_pix_survey = pix_survey#inverted_pix_survey/np.max(inverted_pix_survey)

            levels_ = np.linspace(np.min(inverted_pix_survey),np.percentile(inverted_pix_survey,99),10)

            ax3.contourf(inverted_pix_survey,transform=ax3.get_transform(wcs.WCS(hdr_survey)),levels = levels_,cmap=mp.get_cmap('binary',256))
        except:
            mp.subplot(224)
##
##
##
def PLOT_INDIV_IMG_TESS(clip):

#    disp = lambda x: mp.imshow(x, cmap=mp.cm.YlGnBu_r, origin = "bottom", interpolation = "nearest")#, extent = extent_)
#    disp = lambda x: mp.imshow(x, cmap=parula_map, origin="bottom",interpolation="nearest")

#    itr_mean_img_, oot_mean_img_, diff_mean_img_, cube_itr_mean_, cube_oot_mean_, cube_diff_mean_, transit_number_ = generateImages(clip)
#    itr_mean_img_, oot_mean_img_, diff_mean_img_ = np.asarray(itr_mean_img_), np.asarray(oot_mean_img_), np.asarray(diff_mean_img_)

    out = clip['diffImgCentroids']
    centroids = clip['diffImgCentroids.results']

    ootCol_prf, ootRow_prf = np.mean([centroids[:,0],centroids[:,4]], axis = 0), np.mean([centroids[:,1],centroids[:,5]], axis = 0)
    diffCol_prf, diffRow_prf = centroids[:,2], centroids[:,3]

#    itrCol_psf, itrRow_psf, ootCol_psf, ootRow_psf, diffCol_psf, diffRow_psf = psfCentroids_vbk(clip)
#    diffCol_prf, diffRow_prf = diffCol_psf, diffRow_psf
#    ootCol_prf, ootRow_prf = ootCol_psf, ootRow_psf


    cube = clip['serve.cube']
    time = clip['serve.time']
    period_days = clip['trapFit.period_days']
    epoch_days = clip['trapFit.epoch_bkjd']
    duration_days = clip['trapFit.duration_hrs']/24.

    isnan = np.isnan(time)
    time = time[~isnan]
    cube = cube[~isnan]

    transits = getIngressEgressCadences(time, period_days, epoch_days, duration_days)

    oot_img_ = np.zeros((len(transits), cube.shape[1], cube.shape[2]))
    diff_img_ = np.zeros((len(transits), cube.shape[1], cube.shape[2]))

    transit_flags_ = np.zeros(len(transits), dtype = bool)
    for jj in range(len(transits)):
        key = 'transit-%04i' %(jj)
        if out[key]['errorCode'] < 7:
            cin = transits[jj]
            plot = False
            before, after, diff = generateDiffImg(cube, cin, plot=plot)
            diff_img_[jj,:,:] = diff
            oot_img_[jj,:,:] = 0.5*(before + after)
        else:
            transit_flags_[jj] = True

    transits = transits[~transit_flags_]
    diff_img_, oot_img_ = diff_img_[~transit_flags_], oot_img_[~transit_flags_]
    ootCol_prf, ootRow_prf = ootCol_prf[~transit_flags_], ootRow_prf[~transit_flags_]
    diffCol_prf, diffRow_prf = diffCol_prf[~transit_flags_], diffRow_prf[~transit_flags_]

#    print(oot_img_.shape, oot_img_.shape)
    oot_mean_img_ = np.nanmean(oot_img_, axis = 0)
    diff_mean_img_ = np.nanmean(diff_img_, axis = 0)
    itr_mean_img_ = oot_mean_img_

    ss_ = itr_mean_img_.shape

    extent_ = [0, ss_[1],0, ss_[0]]
#    disp = lambda x: mp.imshow(x, cmap= parula_map, origin = "bottom", interpolation = "nearest", extent = extent_)

#
# PLOT IN-TRANSIT, BEFORE TRANSIT, AND AFTER TRANSIT IMAGES ON A 3x3 GRID	
#
#	subplot_ = 330 + int(ii+1)#10*(int(float(ii)//3.) + 1) + int(ii+1)

    n_columns_ = 3

    transit_number_ = len(transits)

    if transit_number_ > 12:
        transit_number_ = 12

    rows = int(np.ceil(float(transit_number_) / n_columns_))
    gs = gridspec.GridSpec(rows, n_columns_)

    for ii in range(transit_number_):

        mp.subplot(gs[ii])

        plot_ = 'plot_images_'#'plot_flux_'#

        if plot_ == 'plot_flux_':

            number_of_cadences_in_transit_ = transits_[2*ii+1] - transits_[2*ii]

            nn_ = 50
	
            idx_in_transit_ = np.linspace(transits_[2*ii], transits_[2*ii+1], int(number_of_cadences_in_transit_+1))
            idx_in_transit = [int(aa) for aa in idx_in_transit_]

            idx_before_, = np.where(oot_cadence < transits_[2*ii])
            idx_before = oot_cadence[idx_before_[-nn_-1-number_of_cadences_in_transit_:]]

            idx_after_, = np.where(oot_cadence > transits_[2*ii+1])				
            idx_after = oot_cadence[idx_after_[0:number_of_cadences_in_transit_+1+nn_]]
		
            mp.plot(time[idx_in_transit], flux[idx_in_transit], 'r.')
            mp.plot(time[idx_before], flux[idx_before], 'k.')
            mp.plot(time[idx_after], flux[idx_after], 'k.')

        elif plot_ == 'plot_images_':
            kwargs = {'origin':'bottom', 'interpolation':'nearest','cmap':parula_map, 'extent':[0, ss_[1],0, ss_[0]]}
#            disp(cube_diff_mean_[ii,:,:])
            mp.imshow(diff_img_[ii,:,:], **kwargs)
            mp.plot(diffCol_prf[ii], diffRow_prf[ii], 'ro', ms = 6)
            mp.plot(ootCol_prf[ii], ootRow_prf[ii], 'ks', ms = 4)
            mp.colorbar()
#		mp.tight_layout()


def plot_DiffImg_and_Centroids(clip):#cube, centroids, goodCentroidIndices, rollPhase, quality):

    time = clip['serve.time']
    qFlags = clip['serve.flags']

    flux = clip['detrend.flux_frac']
    flags = clip['detrend.flags']

    centroids = clip['diffImg.centroid_timeseries']
    rollPhase = clip['rollPhase.rollPhase']
    period_days = clip['trapFit.period_days']
    epoch_bkjd = clip['trapFit.epoch_bkjd']
    duration_hrs = clip['trapFit.duration_hrs']
    epic = clip['value']

    cube = clip['serve.cube']

    hdr_ = clip['serve.tpfHeader']
    col_zero_, row_zero_ = int(hdr_['1CRV4P']), int(hdr_['2CRV4P'])
    epic_Col, epic_Row = col_zero_ + int(hdr_['1CRPX4']), row_zero_ + int(hdr_['2CRPX4'])

    idx = centroids[:,1] > 0
    cin = centroids[idx, 0]


    ootCol, ootRow = centroids[idx, 1], centroids[idx, 2]
    diffCol, diffRow = centroids[idx, 3], centroids[idx, 4]

    inTransitIndices = kplrfits.markTransitCadences(time, period_days, epoch_bkjd, duration_hrs/24., flags=flags)
    goodCentroidIndices = centroids[ centroids[:,1]>1, 0].asarray().astype(int)

    oot_cadence_, itr_cadence_ = np.where(inTransitIndices == False), np.where(inTransitIndices == True)
    oot_mean_img_, itr_mean_img_ = np.nanmean(cube[oot_cadence_], axis=0), np.nanmean(cube[itr_cadence_], axis=0)

    diff_mean_img_ = oot_mean_img_ - itr_mean_img_

    ss_ = oot_mean_img_.shape

    extent_ = [col_zero_, col_zero_ + ss_[1], row_zero_, row_zero_ + ss_[0]]
    disp = lambda x: mp.imshow(x, cmap=mp.get_cmap('binary', 512), origin = "bottom", interpolation="nearest", extent = extent_)

    n_panels_ = 5*4
    skip_ = 3*15

    for ii in range(n_panels_-1):
        print(ii)
        mp.subplot(5,4,ii+1)

        diff_, oot_, diag = diffimg.constructK2DifferenceImage(cube, goodCentroidIndices[skip_+ii], rollPhase, flags)
        disp(diff_)

        mp.plot(ootCol[skip_+ii], ootRow[skip_+ii], 'c*', ms = 10)
        mp.plot(diffCol[skip_+ii], diffRow[skip_+ii], 'mo', ms = 6)

        mp.title(time[goodCentroidIndices[skip_+ii]], fontsize = 10)
        mp.gca().axes.get_xaxis().set_visible(False)
        mp.gca().axes.get_yaxis().set_visible(False)
        mp.tight_layout()

    mp.subplot(5,4,n_panels_)
    disp(diff_mean_img_)
    mp.scatter(diffCol[skip_+0:skip_+n_panels_], diffRow[skip_+0:skip_+n_panels_], marker='o', c=cin[skip_+0:skip_+n_panels_], s=40, linewidths=0, cmap=mp.cm.RdYlBu)
    cb = mp.colorbar()
    cb.set_label("Cadence")
##
##
##
def plotInTransitAndDiffCentroids(centroids):
    """Plot the row and column of the intranst centroids and the difference
    image centroids"""

    idx = centroids[:,1] > 0
    diffCol = centroids[idx, 'diff_col']
    diffRow = centroids[idx, 'diff_row']

    #itr => in transit
    itrCol = centroids[idx, 'intr_col']
    itrRow = centroids[idx, 'intr_row']

    mp.plot(diffCol, diffRow, 'bo', label="Diff Img Cent")
    mp.plot(itrCol, itrRow, 'rs', label="In-Transit Cent")

    for i in range(len(diffCol)):
        a1 = [diffCol[i], itrCol[i]]
        a2 = [diffRow[i], itrRow[i]]
        mp.plot(a1, a2, '-', color='#AAAAAA', lw=.5)

    mp.legend(loc=0, fancybox=True, shadow=True)
    mp.xlabel(r"Column")
    mp.ylabel(r"Row")
##
##
##
def plotCentroidTimeseries(centroids):
    idx = centroids[:,1] > 0
    rin = centroids[idx, 0]
    diffCol = centroids[idx, 'diff_col']
    diffRow = centroids[idx, 'diff_row']

    #itr => in transit
    itrCol = centroids[idx, 'intr_col']
    itrRow = centroids[idx, 'intr_row']

    diffCol = diffCol - np.mean(diffCol)
    diffRow = diffRow - np.mean(diffRow)

    itrCol = itrCol - np.mean(itrCol)
    itrRow = itrRow - np.mean(itrRow)

#    mp.plot(rin, diffCol, diffRow, 'bo', label="Diff Img Cent")
    ax = mp.subplot(211)
    mp.plot(rin, itrCol,  'rs', label="In-Transit Col")
    mp.plot(rin, diffCol,  'co', label="Diff Img Col")
    mp.ylabel("Deviation from mean Col")

    mp.subplot(212, sharex=ax)
    mp.plot(rin, itrRow,  'rs', label="In-Transit Row")
    mp.plot(rin, diffRow,  'co', label="Diff Img Row")
    mp.ylabel("Deviation from mean Row")

    mp.legend(loc=0)
    mp.xlabel(r"Cadence Number")
##
##
##
def plotCentroidOffsets(centroids):
    """Plot the centroid offsets in column and row

    Inputs:
    ----------
    centroids:
        (Nca) As stored in ``clip['diffImg.centroid_timeseries']``


    Returns:
    ----------
    A string containing the probability the measured centroids are consistent
    with no offset

    Output:
    ------------
    A plot is added to the current figure.
    """
    idx =centroids[:,1] > 0
    cin = centroids[idx, 0]

    ootCol = centroids[idx, 1]
    ootRow = centroids[idx, 2]

    itrCol = centroids[idx, 3]
    itrRow = centroids[idx, 4]

    diffC = -(ootCol - itrCol)
    diffR = -(ootRow - itrRow)

    mp.scatter(diffC, diffR, marker='o', c=cin, s=64, linewidths=0, \
        cmap=mp.cm.RdYlBu)
#    mp.plot(diffC, diffR, 'ko', ms=10)

    mp.axhline(0, color='k', lw=.5)
    mp.axvline(0, color='k', lw=.5)
    mp.plot(0,0, '*', ms=40, color='yellow')

    covar.plotErrorEllipse(diffC, diffR, color='#888888', ms=20)

    mp.xlabel(r"$\Delta$ Column (pixels)")
    mp.ylabel(r"$\Delta$ Row (pixels)")

    probOffset, chiSq = covar.computeProbabilityOfObservedOffset(diffC, diffR)
    titleStr = "Prob. On Target: %.1e: $\chi^2$: %.3f" %(1-probOffset, chiSq)
    cb = mp.colorbar()
    cb.set_label("Cadence")#Time (BKJD)")

#    mp.tight_layout()

    #Ensure some padding around the origin so the symbol is never
    #at edge of the plot
    axl = list(mp.axis())
    axl[0] = min(axl[0], -0.1)
    axl[1] = max(axl[1], +0.1)
    axl[2] = min(axl[2], -0.1)
    axl[3] = max(axl[3], +0.1)
    mp.axis(axl)

    return titleStr
##
##
##
def plotCentroidOffsets_TESS(clip):

    out = clip['diffImgCentroids']
    centroids = clip['diffImgCentroids.results']

    ootCol_prf, ootRow_prf = np.mean([centroids[:,0],centroids[:,4]], axis = 0), np.mean([centroids[:,1],centroids[:,5]], axis = 0)
    diffCol_prf, diffRow_prf = centroids[:,2], centroids[:,3]

#    itrCol_psf, itrRow_psf, ootCol_psf, ootRow_psf, diffCol_psf, diffRow_psf = psfCentroids_vbk(clip)

    diffC, diffR = (ootCol_prf - diffCol_prf), (ootRow_prf - diffRow_prf)

#    itr_mean_img_, oot_mean_img_, diff_mean_img_, itr_mean_cube_, oot_mean_cube_, diff_mean_cube_, transit_number_ = generateImages(clip)
#    itr_mean_img_, oot_mean_img_, diff_mean_img_ = np.asarray(itr_mean_img_), np.asarray(oot_mean_img_), np.asarray(diff_mean_img_)

    hdr_ = clip['serve.tpfHeader']    

    if clip['config.detrendType'] == "eleanor":
    	col_zero_, row_zero_ = int(hdr_['CRPIX1']), int(hdr_['CRPIX2'])
    	epic_Col, epic_Row = col_zero_ + int(hdr_['TPF_W']), row_zero_ + int(hdr_['TPF_H'])

    cube = clip['serve.cube']
    time = clip['serve.time']
    period_days = clip['trapFit.period_days']
    epoch_days = clip['trapFit.epoch_bkjd']
    duration_days = clip['trapFit.duration_hrs']/24.

    isnan = np.isnan(time)
    time = time[~isnan]
    cube = cube[~isnan]

    transits = getIngressEgressCadences(time, period_days, epoch_days, duration_days)

    oot_img_ = np.zeros((len(transits), cube.shape[1], cube.shape[2]))
    diff_img_ = np.zeros((len(transits), cube.shape[1], cube.shape[2]))

    transit_flags_ = np.zeros(len(transits), dtype = bool)
    for jj in range(len(transits)):
        key = 'transit-%04i' %(jj)
        if out[key]['errorCode'] < 7:
            cin = transits[jj]
            plot = False
            before, after, diff = generateDiffImg(cube, cin, plot=plot)
            diff_img_[jj,:,:] = diff
            oot_img_[jj,:,:] = 0.5*(before + after)
        else:
            transit_flags_[jj] = True

    transits = transits[~transit_flags_]
    diff_img_, oot_img_ = diff_img_[~transit_flags_], oot_img_[~transit_flags_]
    ootCol_prf, ootRow_prf = ootCol_prf[~transit_flags_], ootRow_prf[~transit_flags_]
    diffCol_prf, diffRow_prf = diffCol_prf[~transit_flags_], diffRow_prf[~transit_flags_]
    diffC, diffR = diffC[~transit_flags_], diffR[~transit_flags_]

    oot_mean_img_ = np.nanmean(oot_img_, axis = 0)
    diff_mean_img_ = np.nanmean(diff_img_, axis = 0)
    itr_mean_img_ = oot_mean_img_

    ss_ = itr_mean_img_.shape

    extent_ = [0, ss_[1],0, ss_[0]]
    disp = lambda x: mp.imshow(x, cmap= parula_map, origin = "bottom", interpolation = "nearest", extent = extent_)

    disp(oot_mean_img_)

    mp.plot(np.mean(diffCol_prf), np.mean(diffRow_prf), 'r*', ms = 12, mec = 'k', label="Diff")
    mp.plot(diffCol_prf, diffRow_prf, 'r*')#, label="Indiv Diff")
    mp.plot(np.mean(ootCol_prf), np.mean(ootRow_prf), 'mo', label="OOT")
   
    covar.plotErrorEllipse(diffCol_prf, diffRow_prf, color='r', ms=14, marker = '*', mfc = 'r')
    probOffset, chiSq = covar.computeProbabilityOfObservedOffset(diffC, diffR)
    titleStr = "Prob. On Target: %.1e; $\chi^2$: %.3f" %(1-probOffset, chiSq)

    mp.xlabel("Column")
    mp.ylabel("Row")
    mp.legend(loc = 'best', fontsize = 12)
    mp.title('Diff. Image; ' + titleStr, fontsize = 12)

    return titleStr
##
##
##

def multiPanelPlotDiffImgCentroidsDiagnostic(time, flux, flags, rollPhase, \
    inTransitIndices, goodCentroidIndices, qFlags):
    """Plot the flux and roll phase time series.

    A multi panel plot with three panels for the flux time series,
    and three for the roll phase time series. Cadences in transit are marked
    with blue or red squares. Blue indicates a difference image was successfully
    created, red indicatese no difference image.

    Each time series is spread over three panels to better show the details.
    Panels showing common sections of the timeseries are highlighted
    with the same background colour.

    Inputs:
    ----------
    time
        (1d numpy array) X axis for plot

    flux
        (1d numpy array) Flux time series

    flags
        (1d numpy boolean array) Cadences where flags is set are considered
        bad and not plotted.

    rollPhase
        (1d numpy array) Roll phase time series

    inTransitIndices
        (1d numpy boolean array) Which cadences are in transit

    goodCentroidIndices
        (1d numpy boolean array) Which cadences have centroids measured.


    Returns:
    ----------
    **None**

    Output:
    ------------
    A plot is added to the current figure.
    """
    fig = mp.gcf()


    time = np.arange(len(time))

    nPanel = 3
    start = np.min(time[~flags])
    deltaT = np.max(time[~flags]) - start
    deltaT /= float(nPanel)

    colour = ["#FFDDDD", "#DDFFDD", "#DDDDFF"]

    for i in range(nPanel):
        ax = mp.subplot(2*nPanel, 1, 2*i+1, facecolor=colour[i])
        plotTimeseries(time, 1e3*flux, flags, inTransitIndices, \
            goodCentroidIndices, qFlags)
        mp.ylabel("Frac Flux (ppk)")

        mp.subplot(2*nPanel, 1, 2*i+2, sharex=ax, facecolor=colour[i])
        plotTimeseries(time, rollPhase, flags, inTransitIndices, \
            goodCentroidIndices, qFlags)
        mp.ylim(-1.5,1.5)
        mp.ylabel("Roll Phase")

        mp.xlim(start + i*deltaT, start + (i+1)*deltaT)
        mp.xlabel("Cadence")#Time (BKJD)")
##
##
##
#def plotDiffImgCentroidsDiagnostic(time, flux, flags, rollPhase, \
#    inTransitIndices, centroids, qFlags):
#
#    idx =
#
#
##    for i0 in np.where(idx)[0]:
##        i1, i2 = diffimg.getIndicesOfOotImages(rollPhase, i0, flags)
##        indexBefore.append(i1)
##        indexAfter.append(i2)
#
#
#    print centroids[:10]
#
#    mp.clf()
#    ax = mp.subplot(211)
#    plotThrusterFirings(qFlags, time, color='#888888', lw=.4)
#    mp.plot(time[~flags], flux[~flags], 'ko', ms=4, alpha=.8)
#    mp.plot(time[inTransitIndices], flux[inTransitIndices], 'rs')
#    mp.plot(time[goodCentroidIndices], flux[goodCentroidIndices], 'bo')
#
#    mp.subplot(212, sharex=ax)
#    plotThrusterFirings(qFlags, time, color='#888888', lw=.4)
#    mp.plot(time[~flags], rollPhase[~flags], 'ko', ms=4, alpha=.8)
##
##
##
def plotTimeseries(time, y, flags, inTransitIndices, goodCentroidIndices, qFlags):
    """Plot a time series

    Also marks cadences with thruster firings with vertical lines on the plot,
    and the location of cadences in transit, and cadences with measured centroids
    """
    plotThrusterFirings(qFlags, time, color='#888888', lw=.4)
    mp.plot(time[~flags], y[~flags], 'ko', ms=2, alpha=.8)
    mp.plot(time[inTransitIndices], y[inTransitIndices], 'rs')
    mp.plot(time[goodCentroidIndices], y[goodCentroidIndices], 'bo')
##
##
##
def plotThrusterFirings(qualityFlags, xval=None, **kwargs):
    """Mark cadences with thruster firings

    Inputs:
    ----------
    qualityFlags
        (1d numpy array) Quality flags as found in lightcurves and TPF
        files. This is not a boolean array, but an array of flag masks.
        Must include a bit for thruster firings.

    Optional Inputs:
    -----------------
    xval
        (1d array) Values of x axis to plot (eg the time, or cadence number)
        Defeault is 1.. len(qualityFlags)

    All other optional arguments passed to matplotlibs axvline

    Returns:
    ----------
    **None**

    Output:
    ----------
    A series of vertical lines are added to the current plotx
    """

    flags = qualityFlags.astype(np.int32)  #Mneumonic
    if xval is None:
        xval = np.arange(len(flags))

    wh = np.where( flags & kplrfits.SapQuality['DefiniteRollTweak'])[0]

    for w in wh:
        mp.axvline(xval[w], **kwargs)
##
##
##
def plotDiffImg(cube, indexOfCadence, rollPhase, quality):

    if quality.dtype == np.bool:
        raise TypeError("quality should not be a boolean array")

    mp.clf()
    orig = cube[indexOfCadence]
    diff, oot, diag = diffimg.constructK2DifferenceImage(cube, \
        indexOfCadence, \
        rollPhase, quality)
    disp = lambda x: mp.imshow(x, cmap=mp.cm.YlGnBu_r, origin="bottom",interpolation="nearest")
#    disp = lambda x: mp.imshow(x, cmap=parula_map, origin="bottom",interpolation="nearest")

    mp.subplot(221)
    disp(cube[indexOfCadence])
    mp.colorbar()
    mp.title("Cadence")

    mp.subplot(222)
    disp(oot)
    mp.colorbar()
    mp.title("Out of Transit")

    #Quite now if no diff img created
    if np.all(oot == 0):
        return

    mp.subplot(223)
    disp(diff)
    mp.colorbar()
    mp.title("Difference Image")

    mp.subplot(224)
    idx = np.isfinite(diff)
    snr = np.empty_like(diff)
    snr[idx] = diff[idx]/np.sqrt(orig[idx])
    disp(snr)
    mp.colorbar()
    mp.title("Difference SNR")
