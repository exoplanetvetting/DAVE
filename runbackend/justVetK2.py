#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 29 13:44:51 2016

@author: sthomp

Command Line Script to run vetter on the K2 data.
Inputs
File of this info or the info itself
    EpicId
    Campaign
    Period (days)
    epoch (bkjd)
    depth (ppm)
config File 
"""

import dave.pipeline.clipboard as clipboard
import numpy as np
import sys
import os
import getopt as getopt
import dave.pipeline.exporter as ex
import dave.pipeline.multiPagePlot as mpp
import datetime
import dave.pipeline.pipeline as dpp   #You need this
import dave.stellar.readStellarTable as stel

def main():
    """A bare bones main program"""
    

    try:
        opts, args = getopt.getopt(sys.argv[1:], "hf:c:o:1:l:", ["help", "output=","config=","file=","one=","lc="])
    except getopt.GetoptError as err:
        # print help information and exit:
        usage()
        sys.exit()
        
    cfgFile=""
    ephemFile=""
    output=""
    detrendType="pdc"
    data=np.zeros((1,5),dtype=float)  
    print np.shape(data)    
        
    for o, a in opts:
        if o in ("-f","--file"):
            ephemFile = a
            print "Ephemeris File is: %s" % ephemFile
            data=loadEphemFile(ephemFile)
        elif o in ("-h", "--help"):
            usage()
            sys.exit()
        elif o in ("-o", "--output"):
            output = a
            print "Output File is %s\n" % output
        elif o in ("-c", "--config"):
            cfgFile= a
            print "Config File is: %s\n" % cfgFile
        elif o in ("-1", "--one"):
            data[0,:]=np.transpose(np.array(a.split(),dtype=float))
        elif o in  ("-l", "--lc"):
			detrendType=a
        else:
            assert False, "Unhandled option"
            sys.exit()
            
    cfg=loadConfigInput(cfgFile)
    cfg['detrendType']=detrendType
     
    cfg=suppConfiguration(cfg)
    #print cfg 
    
    for i,epic in enumerate(data[:,0]):
        cfg['campaign']=int(data[i,1])
        try:
            dep=data[i,4]/1.0e6
        except:
            dep=.00005
           
        clip=runOneEphem(epic,data[i,2],data[i,3],cfg,duration=3.0,depth=dep) 
        
        print clip.__meta__
        
        if ('exception' not in clip.keys()):
            outfile=runExport(clip,output)
            print 'Created Outputs %s\n\n' % outfile
        else:
            print "No Outputs\n"
            fid=open(output,'a') 
            fid.write("%s %f 0   0   0   0   0   0   0   0 \t-1 -1 -1 -1 NO_Analysis\n" % (clip.value,clip.bls.period))
            fid.close()
            #outfile=runExport(clip,output)
            print clip.exception
            print clip.backtrace
            
 

def usage():
    """Help message
    """
    
    print "justVet -f input ephem file -c config file -o output directory\n"
    print "writes stuff to current directory\n\n"
    print "Format of the input ephem file is\n"
    print "epic campaign period_days epoch depth"
    print "To run just one, use -1 \"epic campaign period epoch depth(ppm)\""
    print "You still need -c and -o"
    print "Use -l or --lc to pick your light curve"
    print "The names of the light curve choices are pdc,everest,sff,agp,varcat"
    print "Default is the PDC light curves."


def loadEphemFile(ephemFile):
    """
    Load a file full of ephemerides
    return data array
    0=epicId
    1=campaign
    2=period
    3=epoch
    4=depth
    """
    
    data=np.loadtxt(ephemFile,dtype=float,comments='#',delimiter=None)
    print "Loaded %s\n" % ephemFile

    return data    
    
    

def loadConfigInput(cfgFile):
    """
    Load a file with information you need for the configuration
    Add into your config file and return.
    """
    cfg={}    
    
    info=np.loadtxt(cfgFile,dtype=str,delimiter=":",comments='#')
    
    for i,key in enumerate(info[:,0]):
        try:
            cfg[key]=float(info[i,1])
        except ValueError:
            cfg[key]=info[i,1]
    
    return cfg
   
   
def suppConfiguration(cfg):
    """Load the default pipeline configuration and adjust as necessary
    """

    #Edit the input configuration with things specific to this task.

    cfg['debug'] = False

#    tasks = """dpp.checkDirExistTask dpp.serveTask dpp.extractLightcurveFromTpfTask
#        dpp.computeCentroidsTask dpp.rollPhaseTask dpp.cotrendSffDataTask
#        dpp.detrendDataTask dpp.trapezoidFitTask dpp.lppMetricTask 
#        dpp.modshiftTask dpp.measureDiffImgCentroidsTask dpp.dispositionTask
#        dpp.saveOnError""".split()  
        
    tasks = """dpp.checkDirExistTask dpp.serveTask dpp.extractLightcurveTask
        dpp.computeCentroidsTask dpp.rollPhaseTask dpp.cotrendDataTask
        dpp.detrendDataTask dpp.trapezoidFitTask dpp.lppMetricTask 
        dpp.modshiftTask dpp.measureDiffImgCentroidsTask dpp.dispositionTask
        dpp.saveClip""".split()
     
    cfg['taskList'] = tasks
    
    searchTaskList = """blsTask trapezoidFitTask modshiftTask
    measureDiffImgCentroidsTask dispositionTask""".split()

    cfg['searchTaskList'] = searchTaskList
    try:
        cfg['timeout_sec'] =  int(cfg['timeout_sec'])
    except:
        cfg['timeout_sec'] = 150
        
        
    return cfg



def runOneEphem(k2id,period,epoch,config,duration=3.5,depth=.0001):
    """
    Run just the vetting and return an output.
    Inputs:
    -------------
    k2id
        (int) Epic id of the target to run on.
    period
        (float) period of the target
    epoch
        (float) Time in days
    config
        (dict) Dictionary of configuration parameters
       
    """
    
    taskList = config['taskList']; 


    clip = clipboard.Clipboard()
    clip['config'] = config
    clip['value'] = k2id
    out = clipboard.Clipboard()
    out['period'] = period
    out['epoch'] = epoch
    out['duration_hrs'] = duration
    out['depth'] = depth
    clip['bls'] = out
            
    

    #Check that all the tasks are properly defined
    for t in taskList:
        f = eval(t)

    #Now run them.
    for t in taskList:
        f = eval(t) 
        clip = f(clip)

    return clip

def runExport(clip,output):
    """
    run the exporters based on the input clip.
    Append the important information to the output File.
    """
    per=np.round(clip.bls.period*10)
    epoch=np.round(clip.bls.epoch)
    
    try:
        clip['config']['stellarPar']=['Mass','Rad','Teff','dis','rho','prov','logg']
        clip=stel.addStellarToClip(clip)
        clip=stel.estimatePlanetProp(clip)
    except:
        print 'No Stellar Values'
        
    outstr,header=ex.createExportString(clip, delimiter=" ", badValue="nan")

    fid=open(output,'a') 
    fid.write("%s\n" % outstr)
    fid.close()    

    tag="%i-%02i-%04i-%s" % (clip.value,per,epoch,clip.config.detrendType)
    outfile="%09i/jvet%s.pdf" % (int(clip.value),tag)

    thedir=str(int(clip.value))
    print thedir
    if ~(os.path.isdir(thedir)):
        try:
            os.makedirs(thedir)
        except:
            print 'Error making directory %s' % clip.value
    
    date=datetime.datetime.now()
    
    if ('disposition' not in clip.keys()):
        clip['disposition'] = 'No Disposition Determined'
        clip.disposition.isCandidate = 0
        clip.disposition.isSignificantEvent = 0
        
    mpp.plot_multipages(outfile, clip, date)

        
    return outfile
    

if __name__ == "__main__":
    main()
