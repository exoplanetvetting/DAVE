import subprocess
import os
import sys
import math
from subprocess import check_output
import time
import numpy



def roboVet(clip):
    """Run the Model-Shift test
    
    Inputs:
    -------------
    clip
        The DAVE clipboard containing all the model-shift results and whatever other tests.
        
    Returns:
    -------------
    A dictionary containing the following keys:
    
    disp
      The disposition of the DOI --- either "candidate" or "false positive"
    not_tran_like
      A 1/0 flag indicating whether or not the DOI is transit-like. 0 means it is transit-like
    sig_sec
      A 1/0 flag indicating whether or not the DOI has a significant secondary.
    comments
      A string containing comments on individual tests
    
    
    Output: 
    ----------
    None
    """
    
    # By default, it is a candidate unless it fails a test.
    disp = 'candidate'
    comments = ''
    
    # Run the not transit-like tests
    out = not_trans_like(clip)
    not_trans_like_flag = out['not_trans_like_flag']
    comments += out['comments']
    
    
    # Run the significant secondary tests
    out = sig_sec(clip)
    sig_sec_flag = out['sig_sec_flag']
    comments += out['comments']
    
    
    # Set false positive if any flag is not 0
    if not_trans_like_flag > 0 or sig_sec_flag > 0:
        disp = 'false positive'
           
        
    return {'disp':disp, 'not_trans_like':not_trans_like_flag, 'sig_sec':sig_sec_flag, 'comments':comments}
  

def not_trans_like(clip):
  
    not_trans_like_flag = 0
    comments = ''
    
    # Check is primary is significant compared to red nosie level
    if clip['modshift']['mod_sig_pri']/clip['modshift']['mod_Fred'] < clip['modshift']['mod_sig_fa1'] and clip['modshift']['mod_sig_pri'] > 0:
        if comments != '':
            comments += '---'
        comments += 'SIG_PRI_OVER_FRED_TOO_LOW'
        not_trans_like_flag = 1

    # Check if primary is significant compared to tertiary
    if clip['modshift']['mod_sig_pri']-clip['modshift']['mod_sig_ter'] < clip['modshift']['mod_sig_fa2'] and clip['modshift']['mod_sig_pri'] > 0 and clip['modshift']['mod_sig_ter'] > 0:
        if comments != '':
            comments += '---'
        comments += 'SIG_PRI_MINUS_SIG_TER_TOO_LOW'
        not_trans_like_flag = 1

    # Check if primary is significant compared to positive
    if clip['modshift']['mod_sig_pri']-clip['modshift']['mod_sig_pos'] < clip['modshift']['mod_sig_fa2'] and clip['modshift']['mod_sig_pri'] > 0 and clip['modshift']['mod_sig_pos'] > 0:
        if comments != '':
            comments += '---'
        comments += 'SIG_PRI_MINUS_SIG_POS_TOO_LOW'
        not_trans_like_flag = 1

    return {'not_trans_like_flag':not_trans_like_flag,'comments':comments}
  
  
def sig_sec(clip):
  
    sig_sec_flag = 0
    comments = ''
    
    # Check if a significant secondary exists in phased light curve from model-shift
    if clip['modshift']['mod_sig_sec'] / clip['modshift']['mod_Fred']    > clip['modshift']['mod_sig_fa1'] and clip['modshift']['mod_sig_sec'] > 0 and \
       clip['modshift']['mod_sig_sec'] - clip['modshift']['mod_sig_ter'] > clip['modshift']['mod_sig_fa2'] and clip['modshift']['mod_sig_ter'] > 0 and \
       clip['modshift']['mod_sig_sec'] - clip['modshift']['mod_sig_pri'] > clip['modshift']['mod_sig_fa2'] and clip['modshift']['mod_sig_pri'] > 0:
        if comments != '':
            comments += '---'
        comments += 'SIG_SEC_IN_MODEL_SHIFT'
        sig_sec_flag = 1
        # Add something here for eclipse from planet refelection one day
        # Next line is to check if it could be detected at twice the orbital period and thus should be a PC.
        # HAVE TO CHECK WITH FERGAL ON VALUES FROM TRAP FIT
        #if abs(0.5 - clip['modshift']['mod_ph_sec'])*clip['trapFit.period_days'] < 0.25*clip['trapFit.duration_hrs']/24.0 and abs(clip['modshift']['mod_sig_pri'] - clip['modshift']['mod_sig_sec']) < clip['modshift']['mod_sig_fa2']
        
    
    # Check Odd/Even from model-shift
    if abs(clip['modshift']['mod_sig_odd'] - clip['modshift']['mod_sig_evn']) > clip['modshift']['mod_sig_fa2']:
        if comments != '':
            comments += '---'
        comments += 'ODD_EVEN_DIFF'
        sig_sec_flag = 1

    
    return {'sig_sec_flag':sig_sec_flag,'comments':comments}
  
  