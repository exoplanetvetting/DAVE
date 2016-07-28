# -*- coding: utf-8 -*-
"""
Created on Thu Jul 28 16:35:57 2016

@author: fergal

$Id$
$URL$
"""


import loadMultipleDetrendings as lmd
import os

def test_smoke():

    epic  = 211816003
    campaign = 5
    dataStorePath = os.path.join(os.environ['HOME'], ",mastio")
    detrendTypes = ["PDC", "Everest", "Agp", "sff"]

    #If this doesn't crash, then all detrendings were loaded.
    return lmd.loadMultipleDetrendings(epic, campaign, dataStorePath,
                                       detrendTypes)