# -*- coding: utf-8 -*-
"""
Created on Tue Feb  26 5:33:35 2019

@author: vkostov

$Id$
$URL$
"""

__version__ = "$Id$"
__URL__ = "$URL$"

import matplotlib.pyplot as mp
import numpy as np
import dave.fileio.kplrfits as kplrfits

def generateImages(clip):

    time = clip['serve.time']
    flux = clip['detrend.flux_frac']
    flags = clip['detrend.flags']

    period_days = clip['trapFit.period_days']
    epoch_bkjd = clip['trapFit.epoch_bkjd']
    duration_hrs = clip['trapFit.duration_hrs']
    epic = clip['value']

    cube = clip['serve.cube']
    hdr_ = clip['serve.tpfHeader']

    inTransitIndices = kplrfits.markTransitCadences(time, period_days, epoch_bkjd, duration_hrs/24., flags=flags)

    oot_cadence_ = np.where((inTransitIndices == False) & (flags == False))
    oot_cadence = np.asarray(oot_cadence_).flatten()

    itr_cadence_ = np.where(inTransitIndices == True)
    itr_cadence = np.asarray(itr_cadence_).flatten()
#
#
# GET IN-TRANSIT, BEFORE TRANSIT, AND AFTER TRANSIT CADENCES ONLY
#
#
    transit_number_ = 1
    transits_ = [itr_cadence[0]]
    tmp_idx_, = np.where(oot_cadence < transits_[0])  
    no_transits_ = oot_cadence[tmp_idx_[-1]]

    for ii in range(1,len(itr_cadence)):
        if itr_cadence[ii] - itr_cadence[ii-1] > 10:
            transit_number_ += 1
            transits_ = np.hstack((transits_ , itr_cadence[ii-1], itr_cadence[ii]))
    transits_ = np.hstack((transits_ , itr_cadence[-1]))

    ss_ = cube.shape
    itr_mean_cube_ = np.zeros((transit_number_, ss_[1], ss_[2]))
    oot_mean_cube_ = np.zeros((transit_number_, ss_[1], ss_[2]))
    diff_mean_cube_ = np.zeros((transit_number_, ss_[1], ss_[2]))

    for ii in range(transit_number_):

        number_of_cadences_in_transit_ = transits_[2*ii+1] - transits_[2*ii]

        idx_in_transit_ = np.linspace(transits_[2*ii], transits_[2*ii+1], int(number_of_cadences_in_transit_+1))
        idx_in_transit = [int(aa) for aa in idx_in_transit_]

        idx_before_, = np.where(oot_cadence < transits_[2*ii])
        idx_before = oot_cadence[idx_before_[-1-number_of_cadences_in_transit_:]]

        idx_after_, = np.where(oot_cadence > transits_[2*ii+1])				
        idx_after = oot_cadence[idx_after_[0:number_of_cadences_in_transit_+1]]

        itr_mean_img_by_transit_ = np.nanmean(cube[idx_in_transit,:,:], axis = 0)
        before_tr_mean_img_by_transit_ = np.nanmean(cube[idx_before,:,:], axis = 0)
        after_tr_mean_img_by_transit_ = np.nanmean(cube[idx_after,:,:], axis = 0)

        oot_mean_img_by_transit_ = 0.5*(before_tr_mean_img_by_transit_ + after_tr_mean_img_by_transit_)

        diff_mean_img_by_transit_ = oot_mean_img_by_transit_ - itr_mean_img_by_transit_

        itr_mean_cube_[ii,:,:] = itr_mean_img_by_transit_
        oot_mean_cube_[ii,:,:] = oot_mean_img_by_transit_
        diff_mean_cube_[ii,:,:] = diff_mean_img_by_transit_

    itr_mean_img_ = np.nanmean(itr_mean_cube_, axis = 0)
    oot_mean_img_ = np.nanmean(oot_mean_cube_, axis = 0)
    diff_mean_img_ = oot_mean_img_ - itr_mean_img_


    return itr_mean_img_, oot_mean_img_, diff_mean_img_, itr_mean_cube_, oot_mean_cube_, diff_mean_cube_, transit_number_
