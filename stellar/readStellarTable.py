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

def addStellarToClip(stellarFile,infoList,clip):
    """
    Load the stellar panda File
    Add a stellar key to the clip and include all the fields in infoList.
    InfoList should have the names as used by the stellar table.
    """
    
    data=pd.read_csv(stellarFile,sep='|',header=3,index_col='#EPIC')
    epic=clip['value']
    new=dict()
    
    for v in infoList:
        new[v]=data.ix[epic][v]
    
    clip['stellar']=new    
    
    return clip
    