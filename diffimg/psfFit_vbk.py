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
from dave.diffimg.generate_images import generateImages

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

    itr_mean_img_, oot_mean_img_, diff_mean_img_, cube_itr_mean_, cube_oot_mean_, cube_diff_mean_, transit_number_ = generateImages(clip)
    itr_mean_img_, oot_mean_img_, diff_mean_img_ = np.asarray(itr_mean_img_), np.asarray(oot_mean_img_), np.asarray(diff_mean_img_)

    itrCol, itrRow, itr_cov = [], [], []
    ootCol, ootRow, oot_cov = [], [], []
    diffCol, diffRow, diff_cov = [], [], []

    cube = clip['serve.cube']
    hdr_ = clip['serve.tpfHeader']

    if clip['config.detrendType'] == "tess_2min":
    	col_zero_, row_zero_ = int(hdr_['1CRV4P']), int(hdr_['2CRV4P'])
    	epic_Col, epic_Row = col_zero_ + int(hdr_['1CRPX4']), row_zero_ + int(hdr_['2CRPX4'])
    
    if clip['config.detrendType'] == "eleanor":
    	col_zero_, row_zero_ = int(hdr_['CRPIX1']), int(hdr_['CRPIX2'])
    	epic_Col, epic_Row = col_zero_ + int(hdr_['TPF_H']), row_zero_ + int(hdr_['TPF_W'])

    cube_back_ = cube
    if clip['config.detrendType'] == "tess_2min":
        cut_ = 3
#	cube = cube[:,epic_Col-col_zero_-cut_:epic_Col-col_zero_+ cut_, epic_Row-row_zero_-cut_:epic_Row-row_zero_+ cut_]
    if clip['config.detrendType'] == "eleanor":
        cut_ = 5
#	cube = cube[:,epic_Col-col_zero_-cut_:epic_Col-col_zero_+ cut_, epic_Row-row_zero_-cut_:epic_Row-row_zero_+ cut_]

#    print cube_itr_mean_[0,:,:].shape
#    print cube_itr_mean_[0,epic_Col-col_zero_-cut_:epic_Col-col_zero_+ cut_, epic_Row-row_zero_-cut_:epic_Row-row_zero_+ cut_].shape
#    print epic_Col, col_zero_
#    print epic_Col-col_zero_, epic_Col-col_zero_- cut_
#    print 

    for ii in range(transit_number_):
		
        itrCol_by_transit_, itrRow_by_transit_, itr_cov_by_transit_ = \
intertial_axis(cube_itr_mean_[ii,epic_Col-col_zero_-cut_:epic_Col-col_zero_+ cut_, epic_Row-row_zero_-cut_:epic_Row-row_zero_+ cut_])
        ootCol_by_transit_, ootRow_by_transit_, oot_cov_by_transit_ = \
intertial_axis(cube_oot_mean_[ii,epic_Col-col_zero_-cut_:epic_Col-col_zero_+ cut_, epic_Row-row_zero_-cut_:epic_Row-row_zero_+ cut_])
        diffCol_by_transit_, diffRow_by_transit_, diff_cov_by_transit_ = \
intertial_axis(cube_diff_mean_[ii,epic_Col-col_zero_-cut_:epic_Col-col_zero_+ cut_, epic_Row-row_zero_-cut_:epic_Row-row_zero_+ cut_])

        itrCol, itrRow = np.hstack((itrCol, itrCol_by_transit_)), np.hstack((itrRow, itrRow_by_transit_))
        ootCol, ootRow = np.hstack((ootCol, ootCol_by_transit_)), np.hstack((ootRow, ootRow_by_transit_))
        diffCol, diffRow = np.hstack((diffCol, diffCol_by_transit_)), np.hstack((diffRow, diffRow_by_transit_))

    return itrCol+(epic_Col-col_zero_-cut_),itrRow+(epic_Row-row_zero_-cut_),ootCol+(epic_Col-col_zero_-cut_),ootRow+(epic_Row-row_zero_-cut_),diffCol+(epic_Col-col_zero_-cut_),diffRow+(epic_Row-row_zero_-cut_)
