import subprocess
import os
import sys
import math
from subprocess import check_output, CalledProcessError
import time
import numpy
import os


def runModShift(time,flux,model,basename,period,epoch):
    """Run the Model-Shift test

    Inputs:
    -------------
    time
        The array of time values days.
    flux
        The array of observed fluxes correspodnding to each time.
    model
        The array of model fluxes corresponding to each time.
    basename
        The basename for the output plot
    period
        The period of the system in days.
    epoch
        The epoch of the system in days.

    Returns:
    -------------
    A dictionary containing the following keys:

    mod_sig_pri
      The significance of the primary event assuming white noise
    mod_sig_sec
      The significance of the secondary event assuming white noise
    mod_sig_ter
      The significance of the tertiary event assuming white noise
    mod_sig_pos
      The significance of the positive event assuming white noise
    mod_sig_odd
      The significance of the primary event utilizing only odd-numbered transits
    mod_sig_evn
      The significance of the primary event utilizing only even-numbered transits
    mod_sig_fa1
      The False Alarm threshold assuming 20,000 objects evaluated
    mod_sig_fa2
      The False Alarm threshold for two events within the phased light curve
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

    Output:
    ----------
    The model-shift plot is also created as a PDF file
    """

#    # Uncomment next 4 lines for testing'mod_sig_fa':mod_sig_fa
#    data = numpy.loadtxt('000757450-01-fulltime-model.dat')
#    time = data[:,0]
#    flux  = data[:,1]
#    model = data[:,2]

    # Write data to a file so it can be read by model-shift compiled C code
    tmpFilename = 'model-shift-in.txt'
    numpy.savetxt(tmpFilename, numpy.c_[time,flux,model])

    # Run modshift, and return the output
    path = getModShiftDir()
    cmd = ["%s/modshift" %(path), 'model-shift-in.txt', basename, \
        str(period), str(epoch)]

    try:
        modshiftcmdout = check_output(cmd)
    except CalledProcessError, e:
        #The called process error message isn't very helpful
        msg = "FAIL: modshift returned error: %s" %(e.output)
        raise IOError(msg)

    # Delete the input text file
    os.remove(tmpFilename)

    # Read the modshift output back in to variables
    info = modshiftcmdout.split()
    mod_sig_pri = float(info[1])
    mod_sig_sec = float(info[2])
    mod_sig_ter = float(info[3])
    mod_sig_pos = float(info[4])
    mod_sig_odd = float(info[5])
    mod_sig_evn = float(info[6])
    mod_sig_fa1 = float(info[7])
    mod_sig_fa2 = float(info[8])
    mod_Fred    = float(info[9])
    mod_ph_pri  = float(info[10])
    mod_ph_sec  = float(info[11])
    mod_ph_ter  = float(info[12])
    mod_ph_pos  = float(info[13])
    mod_secdepth = float(info[14])
    mod_secdeptherr = float(info[15])


    return {'mod_sig_pri':mod_sig_pri, 'mod_sig_sec':mod_sig_sec, 'mod_sig_ter':mod_sig_ter, 'mod_sig_pos':mod_sig_pos, 'mod_sig_odd':mod_sig_odd, 'mod_sig_evn':mod_sig_evn, 'mod_sig_fa1':mod_sig_fa1, 'mod_sig_fa2':mod_sig_fa2, 'mod_Fred':mod_Fred, 'mod_ph_pri':mod_ph_pri, 'mod_ph_sec':mod_ph_sec, 'mod_ph_ter':mod_ph_ter, 'mod_ph_pos':mod_ph_pos, 'mod_secdepth':mod_secdepth, 'mod_secdeptherr':mod_secdeptherr}



def getModShiftDir():
    """Get the path where Mod shift stores its .m files"""
    pathSep = "/"
    path = os.path.realpath(__file__)
    path = pathSep.join(path.split(pathSep)[:-1])
    return path
