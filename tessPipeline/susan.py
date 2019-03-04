# -*- coding: utf-8 -*-
# Copyright 2017-2018 Orbital Insight Inc., all rights reserved.
# Contains confidential and trade secret information.
# Government Users: Commercial Computer Software - Use governed by
# terms of Orbital Insight commercial license agreement.

"""
Created on Tue Nov 27 21:25:04 2018

astropy
scikit-learn
lpproj (pip)

@author: fergal
"""

#from __future__ import print_function
from __future__ import division

from pdb import set_trace as debug
import pandas as p
import numpy as np
import sys
import os
import getopt as getopt
import dave.pipeline.multiPagePlot as mpp
import dave.tessPipeline.tessfunc as tessfunc
import dave.pipeline.clipboard as clipboard
from dave.pipeline.task import task
import dave.tessPipeline.tessPipeline as tessPipeline
import dave.vetting.RoboVet as RoboVet
import dave.pipeline.exporter as ex
import datetime
import dave.stellar.readStellarTable as stel


def createConfig(sector, tic, planetNum, debugMode=True):
    cfg = dict()
    cfg['debug'] = debugMode
    cfg['sector'] = sector
    cfg['campaign'] = sector
    cfg['tic'] = tic
    cfg['planetNum'] = planetNum
    cfg['value'] = tic
    cfg['detrendType'] = "tess"

    #TODO This shouldn't be hardcoded, but passed as a parameter
    cfg['dvtLocalPath'] = "/Users/vkostov/Desktop/Ideas_etc/DAVE_test/TESSting/tmp/"
    
    #TODO Need modshift paths
    cfg['lppMapFile'] = "/Users/vkostov/Desktop/Ideas_etc/DAVE_test/TESSting/LPP_map/combMapDR25AugustMapDV_6574.mat"

    cfg['modshiftBasename'] = "/Users/vkostov/Desktop/Ideas_etc/DAVE_test/TESSting/tmp"       
    cfg['onepageBasename'] = "/Users/vkostov/Desktop/Ideas_etc/DAVE_test/TESSting/justVet/"

    cfg['taskList'] = ['serveTask', 'trapezoidFitTask','modshiftTask', 'sweetTask', 'lppMetricTask']
    #'computeCentroidsTask', 'rollPhaseTask', 
    clip = clipboard.Clipboard(cfg)
    
    return clip
    

def outputInfo(clip):
    """Createa  text string to output to a file
    """     
    header="tic,planetNum"
    
    text = "%i" % (clip.config.tic)
    text = "%s,%i" % (text, clip.config.planetNum)
    
    for k in clip.serve.param.keys():
        text = "%s,%f" % (text, clip['serve']['param'][k])
        header = "%s,%s" % (header, k)

    text = "%s,%f" % (text, clip.lpp.TLpp)
    header = "%s,Tlpp" % (header) 
    
    if "modshiftTask" in clip.config.taskList:
        for k in clip.modshift.keys():    
            text = "%s,%f" % (text, clip['modshift'][k])
            header = "%s,%s" % (header, k)
    
    for k in clip.robovet.keys():
        text = "%s, %s" % (text, str(clip['robovet'][k]))
        header = "%s,rv_%s" % (header, k)
    
    return text,header


def runOneDv(sector,tic,planetNum,debugMode=True):
    
    cfg = createConfig(sector,tic,planetNum,debugMode=True)
    
    clip = tessPipeline.runOne(cfg,returnClip=True)
    
    out = clipboard.Clipboard(isSignificantEvent=True, isCandidate=True, reasonForFail="None")
    # This is an attempt at quickly vetting the signals.
    # This should be its own task.
    try:
        rv = RoboVet.roboVet(clip.modshift)

	out['fluxVet'] = rv
	clip['disposition'] = out
	clip['disposition.isSignificantEvent']
    	clip['disposition.isCandidate']
    	clip['disposition.reasonForFail']

    except:
        rv={}
        rv['disp'] = 'candidate'
        rv['comments'] = 'NO_MODSHIFT'
        rv['not_trans_like'] = 0
        rv['sig_sec'] = 0
    
    lpp_th = 4  #Threshold for LPP
    
    if clip.lpp.TLpp > lpp_th:
        rv['disp'] = 'false positive'
        rv['comments']= rv['comments'] + "-LPP_TOO_HIGH"
        rv['not_trans_like']=1
    
    sweet_th = 3.5
    if (clip.sweet.amp[0,-1] > sweet_th) | \
       (clip.sweet.amp[1,-1] >sweet_th) | \
       (clip.sweet.amp[2,-1] > sweet_th):
           rv['disp']='false positive'
           rv['comments']=rv['comments'] + "-SWEET_FAIL"
           rv['not_trans_like']=1
    
    clip['robovet']=rv
    
    return clip
    
def runAllTces(tceFile,sector,outfile):
    """
    Run for all TCEs
    """
    df=p.read_csv(tceFile,comment='#')
    
    for index,row in df[0:2].iterrows():
        clip = runOneDv(sector,row.ticid,row.planetNumber)
        text,header = outputInfo(clip)
        
        #Write out decision
        with open(outfile,"a") as fp:
        
            if index == 0:
                fp.write(header+"\n")
                
            fp.write(text + "\n")

	print row.ticid, row.planetNumber, row.orbitalPeriodDays, row. transitEpochBtjd, "... DONE!"


    return clip

def runExport(clip,output):
    """
    run the exporters based on the input clip.
    Append the important information to the output File.
    """
    per=np.round(clip['serve.param.orbitalPeriod_days']*10)
    epoch=np.round(clip['serve.param.epoch_btjd'])
    basedir=clip.config['onepageBasename']
#    try:
#        clip['config']['stellarPar']=['Mass','Rad','Teff','dis','rho','prov','logg']
#        clip=stel.addStellarToClip(clip)
#        clip=stel.estimatePlanetProp(clip)
#    except:
#        print 'No Stellar Values'
        
    outstr,header=ex.createExportString(clip, delimiter=" ", badValue="nan")

    fid=open(output,'a') 
    #fid.write("%s\n" % header)
    fid.write("%s\n" % outstr)
    fid.close()    

    tag="%i-%02i-%04i-%s" % (clip.config.value,per,epoch,clip.config.detrendType)
    outfile="%s/%09i/jvet%s" % (basedir,int(clip.config.value),tag)
    
    thedir=basedir + str(int(clip.config.value))
    try:
        os.mkdir(thedir)
    except OSError:
        donothing = -999.
#        print "Cannot create directory " + thedir
    #print thedir
    
    date=datetime.datetime.now()

    if ('disposition' not in clip.keys()):
        clip['disposition'] = 'No Disposition Determined'
        clip.disposition.isCandidate = 0
        clip.disposition.isSignificantEvent = 0
       
    clip['value'] = clip.config.value

#    mpp.plot_multipages(outfile, clip, date)

        
    return outfile
