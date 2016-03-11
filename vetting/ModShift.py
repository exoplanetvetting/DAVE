from subprocess import check_output, CalledProcessError
import tempfile
import numpy
import os


def runModShift(time,flux,model,plotname,objectname,period,epoch):
    """Run the Model-Shift test

    Inputs:
    -------------
    time
        The array of time values days.
    flux
        The array of observed fluxes correspodnding to each time.
    model
        The array of model fluxes corresponding to each time.
    plotname
        The name for the output plot
    objectname
        The name of the object, to be displayed in the plot title
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
    mod_sig_oe
      The significance of the odd-even metric
    mod_dmm
      The ratio of the individual depths's median and mean values.
    mod_shape
      The shape metric.
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


    timeout_sec = 10
    # Write data to a tempoary file so it can be read by model-shift
    # compiled C code. mkstemp ensures the file is written to a random
    # location so the code can be run in parallel.
    tmpFilename = tempfile.mkstemp(prefix="modshift-%s" %(objectname))[1]
    numpy.savetxt(tmpFilename, numpy.c_[time,flux,model])

    # Run modshift, and return the output
    #Python 2's subprocess module does not easily support timeouts, so
    #instead we use the shell's version
    #Insert rant here asking why subprocess doesn't have a timeout when it's
    #the complicated module that was supposed to communication better.
    path = getModShiftDir()
    cmd = ["timeout", "%i" %(timeout_sec),  "%s/modshift" %(path), \
       tmpFilename, plotname, objectname, str(period), str(epoch)]

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
    mod_sig_oe = float(info[5])
    mod_dmm = float(info[6])
    mod_shape = float(info[7])
    mod_sig_fa1 = float(info[8])
    mod_sig_fa2 = float(info[9])
    mod_Fred    = float(info[10])
    mod_ph_pri  = float(info[11])
    mod_ph_sec  = float(info[12])
    mod_ph_ter  = float(info[13])
    mod_ph_pos  = float(info[14])
    mod_secdepth = float(info[15])
    mod_secdeptherr = float(info[16])


    return {'mod_sig_pri':mod_sig_pri, 'mod_sig_sec':mod_sig_sec, 'mod_sig_ter':mod_sig_ter, 'mod_sig_pos':mod_sig_pos, 'mod_sig_oe':mod_sig_oe, 'mod_dmm':mod_dmm, 'mod_shape':mod_shape, 'mod_sig_fa1':mod_sig_fa1, 'mod_sig_fa2':mod_sig_fa2, 'mod_Fred':mod_Fred, 'mod_ph_pri':mod_ph_pri, 'mod_ph_sec':mod_ph_sec, 'mod_ph_ter':mod_ph_ter, 'mod_ph_pos':mod_ph_pos, 'mod_secdepth':mod_secdepth, 'mod_secdeptherr':mod_secdeptherr}



def getModShiftDir():
    """Get the path where Mod shift stores its .m files"""
    pathSep = "/"
    path = os.path.realpath(__file__)
    path = pathSep.join(path.split(pathSep)[:-1])
    return path
