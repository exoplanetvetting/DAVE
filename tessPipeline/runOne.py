#!/usr/bin/env python
# -*- coding: utf-8 -*-

from dave.tessPipeline import vet_tess_ as pp_tess

detrendType = "eleanor"#"tess"#

#clip = pp_tess.runOneDv("tess", 5, 307210830, 2, 2.2533, 1438.26950, 787, 1.7)
clip = pp_tess.runOneDv(detrendType, 1,271893367,1,5.8707,1326.2738,5338.0,2.4)
pp = pp_tess


#if detrendType == "tess":
#    clip = pp_tess.runOneDv(detrendType, 2, 307210830, 1, 3.690613, 1356.2038, 1863, 1.27)
#    clip = pp_tess.runOneDv(detrendType, 5, 307210830, 1, 3.69046, 1441.0885, 1863, 1.27)
#    clip = pp_tess.runOneDv(detrendType, 2, 307210830, 2, 7.451113, 1355.2864, 1641, 0.8)
#    clip = pp_tess.runOneDv(detrendType, 5, 307210830, 2, 7.451113, 1444.6963, 1641, 0.8)
#    clip = pp_tess.runOneDv(detrendType, 2, 307210830, 3, 2.2530, 1354.90621, 571, 1.)
#    clip = pp_tess.runOneDv(detrendType, 5, 307210830, 2, 2.2533, 1438.26950, 787, 1.)

#    pp = pp_tess
#else:
#    clip = pp_eleanor.runOneDv(detrendType,1,142018015,1,2.7355,1.66,1.34,2.4)
#    clip = pp_eleanor.runOneDv(detrendType,1,149038231,1,11.256,11.235,1.13,2.4)
#    clip = pp_eleanor.runOneDv(detrendType,1,150513899,1,3.7734,0.245,10.94,2.4)
#    clip = pp_eleanor.runOneDv(detrendType,1,31961007,1,4.2554,1.745,107.78,2.4)
#    clip = pp_eleanor.runOneDv(detrendType,1,38462522,1,2.7154,1.865,2.88,2.4)
#    clip = pp_eleanor.runOneDv(detrendType,1,38827240,1,2.3747,1.43,14.76,2.4)
#    clip = pp_eleanor.runOneDv(detrendType,1,141811446,1,3.5715,0.8,11.77,2.4)
#    clip = pp_eleanor.runOneDv(detrendType,1,177162886,1,12.843,1.16,6.72,2.4)
#    clip = pp_eleanor.runOneDv(detrendType,1,271893367,1,5.873,0.95,5.43,2.4)
#    clip = pp_eleanor.runOneDv(detrendType,1,388128415,1,7.3222,1.67,1.98,2.4)
#    clip = pp_eleanor.runOneDv(detrendType,1,308544271,1,2.4959,0.545,19.18,2.4)
#    clip = pp_eleanor.runOneDv(detrendType,1,404851003,1,3.0253,2.615,0.98,2.4)
#    clip = pp_eleanor.runOneDv(detrendType,1,149992242,1,2.6658,0.85,1.94,2.4)
#    clip = pp_eleanor.runOneDv(detrendType,1,177283727,1,3.956,1.735,5.39,2.4)
#XXXXX
#    clip = pp_eleanor.runOneDv(detrendType,1,176871438, 1, 15.0892, 11.465, 8060, 2.4)
#    clip = pp_eleanor.runOneDv(detrendType,1,141768070, 1, 2.786, 1327.19, 11900, 0.5)
#    clip = pp_eleanor.runOneDv(detrendType,1,177283727, 1, 3.957, 1327.0635, 6390, 4.01)
#    clip = pp_eleanor.runOneDv(detrendType,1,30848598, 1, 0.754, 1325.68, 540, 0.5)
#    clip = pp_eleanor.runOneDv(detrendType,1,31415158, 1, 0.790, 1326.05, 190, 0.5)
#    clip = pp_eleanor.runOneDv(detrendType,1,140691463, 1, 2.084, 1326.55, 5100, 0.3)
#    clip = pp_eleanor.runOneDv(detrendType,1,279322914, 1, 18.866, 1337.66, 10760, 6.42)
#    clip = pp_eleanor.runOneDv(detrendType,1,260304800, 1, 3.441, 1327.4983, 1470, 0.5)
#    clip = pp_eleanor.runOneDv(detrendType,1,30848598, 1, 0.754, 1325.6466, 840, 2.39)
#    clip = pp_eleanor.runOneDv(detrendType,1,31415158, 1, 0.790, 1326.0621, 990, 1.94)
#    clip = pp_eleanor.runOneDv(detrendType,1,38462522, 1, 1.357, 1325.8451, 3980, 2.91)
#    clip = pp_eleanor.runOneDv(detrendType,1,167343388, 1, 2.185, 1325.4504, 2070, 2.36)
#    clip = pp_eleanor.runOneDv(detrendType,1,260304800, 1, 3.441, 1327.4983, 3470, 2.39)
#    clip = pp_eleanor.runOneDv(detrendType,1,279322914, 1, 18.866, 1337.6560, 10760, 6.42)
#    clip = pp_eleanor.runOneDv(detrendType,1,141768070, 1, 2.786, 1327.2024, 11900, 0.5)
#XXXXX
#    clip = pp_eleanor.runOneDv(detrendType,1,372909935, 1, 3.611, 1328.5028, 8290, 1.75)
#    clip = pp_eleanor.runOneDv(detrendType,1,38462522,1,1.3559,1325.8642,2198.0,1.5)
#    clip = pp_eleanor.runOneDv(detrendType,1,141811446,1,1.7851,1326.1444,11152.0,1.5)
#    clip = pp_eleanor.runOneDv(detrendType,1,149601448,1,0.8855,1325.853,5391.0,1.5)
#    clip = pp_eleanor.runOneDv(detrendType,1,149990841,1,4.4577,1328.0233,3798.0,1.5)
#    clip = pp_eleanor.runOneDv(detrendType,1,149992242,1,1.3328,1326.1734,2011,2.4)
#    clip = pp_eleanor.runOneDv(detrendType,1,150101799,1,0.7882,1325.4432,3657.0,2.4)
#    clip = pp_eleanor.runOneDv(detrendType,1,150513899,1,1.8857,1325.5783,10620.0,2.4)
#    clip = pp_eleanor.runOneDv(detrendType,1,167366701,1,2.228,1326.3283,11815.0,2.4)
#    clip = pp_eleanor.runOneDv(detrendType,1,177115195,1,14.441,1339.7136,9307.0,2.4)
#    clip = pp_eleanor.runOneDv(detrendType,1,271893367,1,5.8707,1326.2738,5338.0,2.4)
#    clip = pp_eleanor.runOneDv(detrendType,1,306576384,1,3.2151,1327.8433,1926.0,2.4)
#    clip = pp_eleanor.runOneDv(detrendType,1,177283727,1,3.9561,1327.0589,5380.0,2.4)
#    clip = pp_eleanor.runOneDv(detrendType,1,260536927,1,5.7248,1330.1427,7623.0,2.4)
#    clip = pp_eleanor.runOneDv(detrendType,1,300244163,1,2.1799,1326.4383,14260.0,2.4)
#    clip = pp_eleanor.runOneDv(detrendType,1,364395234,1,1.3755,1326.6026,4400.0,2.4)
#    clip = pp_eleanor.runOneDv(detrendType,1,364424328,1,0.7292,1325.8146,15113.0,2.4)
#    clip = pp_eleanor.runOneDv(detrendType,1,375037388,1,0.613,1325.3301,4045.0,2.4)
#    pp = pp_eleanor

outfile_ = 'tmp.txt'
aa = pp_tess.runExport(clip,outfile_)

#clip = susan.runOneDv(1, 38462193, 1)
#clip = susan.createConfig(2, 307210830, 1)
#aa = susan.outputInfo(clip)

#https://mast.stsci.edu/api/v0.1/Download/file/?uri=mast:TESS/product/tess2018206045859-s0001-0000000300903537-0120-s_lc.fits