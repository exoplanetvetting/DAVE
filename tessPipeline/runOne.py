#!/usr/bin/env python
# -*- coding: utf-8 -*-

from dave.tessPipeline import vet_tess_ as pp_tess

detrendType = "eleanor"#"tess"#

#clip = pp_tess.runOneDv(detrendType, 5, 307210830, 2, 2.2533, 1438.26950, 787, 1.7)
clip = pp_tess.runOneDv(detrendType, 1,271893367,1,5.8707,1326.2738,5338.0,2.4)
pp = pp_tess

outfile_ = 'tmp.txt'
aa = pp_tess.runExport(clip,outfile_)
