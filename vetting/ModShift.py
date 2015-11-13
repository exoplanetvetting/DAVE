import subprocess
import os
import sys
import math
from subprocess import check_output
import time
import numpy


def runModShift(phase,flux,model,basename,period):
    """Run the Model-Shift test
    Inputs:
    -------------
    phase
        The array of phases in days, ranging from -0.25*period to 0.75*period.
    flux
        The array of observed fluxes
    model
        The array of model fluxes
    basename
        The basename for the output plot
    period
        The period of the system in days
        
    Returns:
    
    mod_sig_pri
      The significance of the primary event assuming white noise
    mod_sig_sec
      The significance of the secondary event assuming white noise
    mod_sig_ter
      The significance of the tertiary event assuming white noise
    mod_sig_pos
      The significance of the positive event assuming white noise
    mod_sig_fa
      The False Alarm threshold assuming 20,000 objects evaluated
    mod_Fred
      The ratio of the red noise to the white noise in the phased light curve at the transit timescale
    mod_ph_pri
      The phase of the primary event
    mod_ph_sec
      The phase of the secondary event
    mod_ph_sec
      The phase of the tertiary event
    mod_ph_pos
      The phase of the primary event
    mod_secdepth
      The depth of the secondary event
    mod_secdeptherr
      The error in the depth of the secondary event
   
    -------------
    **None**
    Output: 
    ----------
    The model-shift plot is also created as a PDF file
    """
    
    # Uncomment for testing
    data = numpy.loadtxt('000757450-01.mod')
    phase = data[:,0]
    flux  = data[:,1]
    model = data[:,2]
    
    # Write data to a file so it can be read by model-shift compiled C code
    numpy.savetxt('model-shift-in.txt', numpy.c_[phase,flux,model])
    
    # Run modshift, and return the output
    modshiftcmdout = check_output(["./modshift",'model-shift-in.txt',basename,str(period)])
    
    # Delete the input text file
    os.remove('model-shift-in.txt')
    
    # Read the modshift output back in to variables
    info = modshiftcmdout.split()
    mod_sig_pri = float(info[1])
    mod_sig_sec = float(info[2])
    mod_sig_ter = float(info[3])
    mod_sig_pos = float(info[4])
    mod_sig_fa = float(info[5])
    mod_Fred = float(info[6])
    mod_ph_pri = float(info[7])
    mod_ph_sec = float(info[8])
    mod_ph_ter = float(info[9])
    mod_ph_pos = float(info[10])
    mod_secdepth = float(info[11])
    mod_secdeptherr = float(info[12])
    
    
    return {'mod_sig_pri':mod_sig_pri, 'mod_sig_sec':mod_sig_sec, 'mod_sig_ter':mod_sig_ter, 'mod_sig_pos':mod_sig_pos, 'mod_sig_fa':mod_sig_fa, 'mod_Fred':mod_Fred, 'mod_ph_pri':mod_ph_pri, 'mod_ph_sec':mod_ph_sec, 'mod_ph_ter':mod_ph_ter, 'mod_ph_pos':mod_ph_pos, 'mod_secdepth':mod_secdepth, 'mod_secdeptherr':mod_secdeptherr}


