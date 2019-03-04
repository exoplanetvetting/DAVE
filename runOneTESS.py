#!/usr/bin/env python
# -*- coding: utf-8 -*-

from dave.tessPipeline import vet_tess_

detrendType_ = "eleanor"#"tess"#

# Input is: Sector, TIC ID, Planet Number, Period, BTJD, Transit Depth [ppm], Transit Duration [hours]
clip = vet_tess_.runOneDv(detrendType_, 307210830, 1, 3.690613, 1356.2038, 1863, 1.27)
#clip = vet_tess_.runOneDv(detrendType_, 1, 271893367,1,5.8707,1326.2738,5338.0,2.4)

outfile_ = 'tmp.txt'
export_to_file_ = vet_tess_.runExport(clip,outfile_)
