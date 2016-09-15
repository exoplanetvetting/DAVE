# -*- coding: utf-8 -*-
"""
Created on Fri Feb 19 11:27:02 2016

@author: smullall
"""

import numpy as np
import pickle as pk
import dave.fileio.nca as nca
import pandas as pd
import dave.pipeline.clipboard as clipboard

def createStellarNca(stellarFile):
    """
    Read in the stellar file and return a dictionary
    """
    epics=np.loadtxt(stellarFile,dtype=str,delimiter='|',usecols=[0])
    
    data=np.loadtxt(stellarFile,dtype=float,delimiter='|',usecols=np.arange(1,25))
    col='Teff|sTeff+|sTeff-|logg|slogg+|slogg-|FeH|sFeH+|sFeH-|Rad|sRad+|sRad-|Mass|sMass+|sMass-|rho|srho+|srho-|dis|dis+|dis-|ebv|ebv+|ebv-'.split('|');
    
    nameDict=dict()
    nameDict[0]=list(epics)
    nameDict[1]=col
    
    d=nca.Nca(data,nameDict)
    
    return d

def createStellarPanda(stellarFile):    
    """
    Read in the stellar file and return a dictionary -- spaces were removed from file
    data.ix[201121245]['Teff'] is example how to get vlaues out
    """
  
    data=pd.read_csv(stellarFile,sep='|',header=3,index_col='#EPIC')

    return data

def addStellarToClip(clip):
    """
    Load the stellar panda File
    Add a stellar key to the clip and include all the fields in infoList.
    InfoList should have the names as used by the stellar table.
    """
    
    infoList=clip.config.stellarPar
    stellarFile=clip.config.stellarFile

    data=pd.read_csv(stellarFile,sep='|',header=3,index_col='#EPIC')
    epic=int(clip['value'])
    new=dict()
    
    
    for v in infoList:
        try:
            new[v]=data.ix[epic][v]
            
        except KeyError,e:
            new[v]=data.ix[0][v]
            
        clip['stellar']=new    
    
    return clip
    
def estimatePlanetProp(clip):
    """
    Given the information in the clip, depth adn stellar info
    Return a planet dictionary with useful information
    like planet radii and insolation flux.
    Must run addStellarToClip with the following added
    Mass
    Rad
    Teff
    """
    
    srad=clip.stellar.Rad
    stemp=clip.stellar.Teff
    smass=clip.stellar.Mass
    depth=clip.trapFit.depth_frac
    period_sec=clip.trapFit.period_days*24*3600;
    G_si=6.67408e-11
    AU=1.4960e11
    solarR_km=6.957e5
    solarM_kg=1.98855e30
    earthR_km=6371.0
    solarR_earth=solarR_km/earthR_km

    new=dict()
    
    prad_km=solarR_km*srad*(depth**(0.5))

    new['rad_earth']=prad_km/earthR_km;
    
    #Calculate semi-major axis
    acubed=G_si*(period_sec**(2))*solarM_kg*smass/(4*(np.pi)**2)
    a=acubed**(0.33333)

    new['sma_m']=a
    new['sma_au']=a/AU
    
    clip['planet']=new
    
    return clip
    
    
    
    
    
    
    
    