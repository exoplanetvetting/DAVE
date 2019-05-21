#!/usr/bin/env python
# -*- coding: utf-8 -*-

from dave.tessPipeline import vet_tess_

detrendType_ = "tess_2min"#"eleanor"#
sector_ = 2
tic_ = 307210830
planet_num_ = 1
Period_ = 3.690613 # Days
bjd_ = 1356.2038
tdep_ = 1863 # transit depth in ppm
tdur_ = 1.27 # transit duration in hours

# Input is: Sector, TIC ID, Planet Number, Period, BTJD, Transit Depth [ppm], Transit Duration [hours]
clip = vet_tess_.runOneDv(detrendType_, sector_, tic_, planet_num_, Period_, bjd_, tdep_, tdur_)
#clip = vet_tess_.runOneDv(detrendType_, 1, 271893367,1,5.8707,1326.2738,5338.0,2.4)

outfile_ = 'tmp.txt'
export_to_file_ = vet_tess_.runExport(clip,outfile_)
