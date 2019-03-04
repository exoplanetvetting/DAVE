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

def raw_moment(data, iord, jord):
    nrows, ncols = data.shape
    y, x = np.mgrid[:nrows, :ncols]
    data = data * x**iord * y**jord
    return data.sum()

def intertial_axis(data):
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

def make_lines(eigvals, eigvecs, mean, i):
    std = np.sqrt(eigvals[i])
    vec = 3 * std * eigvecs[:,i] / np.hypot(*eigvecs[:,i])
    x, y = np.vstack((mean-vec, mean, mean+vec)).T
    return x, y

def psfCentroids_vbk(clip):

    itrCol, itrRow, itr_cov = [], [], []
    ootCol, ootRow, oot_cov = [], [], []
    diffCol, diffRow, diff_cov = [], [], []

    time = clip['serve.time']
    flux = clip['detrend.flux_frac']
    flags = clip['detrend.flags']

    period_days = clip['trapFit.period_days']
    epoch_bkjd = clip['trapFit.epoch_bkjd']
    duration_hrs = clip['trapFit.duration_hrs']
    epic = clip['value']

    cube = clip['serve.cube']
    hdr_ = clip['serve.tpfHeader']

    if clip['config.detrendType'] == "tess":
    	col_zero_, row_zero_ = int(hdr_['1CRV4P']), int(hdr_['2CRV4P'])
    	epic_Col, epic_Row = col_zero_ + int(hdr_['1CRPX4']), row_zero_ + int(hdr_['2CRPX4'])
    
    if clip['config.detrendType'] == "eleanor":
    	col_zero_, row_zero_ = int(hdr_['CRPIX1']), int(hdr_['CRPIX2'])
    	epic_Col, epic_Row = col_zero_ + int(hdr_['TPF_H']), row_zero_ + int(hdr_['TPF_W'])

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

    cube_back_ = cube
    if clip['config.detrendType'] == "tess":
	cut_ = 3
	cube = cube[:,epic_Col-col_zero_-cut_:epic_Col-col_zero_+ cut_, epic_Row-row_zero_-cut_:epic_Row-row_zero_+ cut_]
    if clip['config.detrendType'] == "eleanor":
	cut_ = 5
	cube = cube[:,epic_Col-col_zero_-cut_:epic_Col-col_zero_+ cut_, epic_Row-row_zero_-cut_:epic_Row-row_zero_+ cut_]

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
		
	itrCol_by_transit_, itrRow_by_transit_, itr_cov_by_transit_ = intertial_axis(itr_mean_img_by_transit_)
	ootCol_by_transit_, ootRow_by_transit_, oot_cov_by_transit_ = intertial_axis(oot_mean_img_by_transit_)
	diffCol_by_transit_, diffRow_by_transit_, diff_cov_by_transit_ = intertial_axis(diff_mean_img_by_transit_)

	itrCol, itrRow = np.hstack((itrCol, itrCol_by_transit_)), np.hstack((itrRow, itrRow_by_transit_))
	ootCol, ootRow = np.hstack((ootCol, ootCol_by_transit_)), np.hstack((ootRow, ootRow_by_transit_))
	diffCol, diffRow = np.hstack((diffCol, diffCol_by_transit_)), np.hstack((diffRow, diffRow_by_transit_))


    return cut_+itrCol, cut_+itrRow, cut_+ootCol, cut_+ootRow, cut_+diffCol, cut_+diffRow