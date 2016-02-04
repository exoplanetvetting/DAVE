""" Module to perform a trapezoid model fit to flux time seres data
    Author: Christopher J Burke
"""
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt

def phaseData(t, per, to):
    """Phase the data at period per and centered at to
       INPUT:
         t - time of data
         per - period to phase time data period and t should
               be in same units
         to - epoch of phase zero
       OUTPUT:
         phi - data phased running from -0.5<phi<=0.5
     """
    phi = np.mod(t - to, per) / per
    phi = np.where(phi > 0.5, phi - 1.0, phi)
    return phi

class trp_parameters:
    """Storage class for the parameters of the trapezoid fit algorithms
       CONTENTS:
         samplen - [Int] subsampling of LC model data
                   ***MUST BE ODD***  No checking this
         likehoodmoddisplay - [Int] If debugLevel > =3 display likelihood call
                              model and residual every iteration mod of
                              this parameter
         cadlen - [Days] Cadence duration
         fitregion - [float] Factor of duration around midpoint to actually
                     fit to data.
    """
    def __init__(self):
        self.samplen = 15
        self.likehoodmoddisplay = 10
        self.cadlen = 29.424/60.0/24.0 #Kepler cadence
        self.fitregion = 4.0
        self.debugLevel = 4

    def __str__(self):
        for k in self.__dict__:
            print k, self.__dict__[k]
        return ''

class trp_originalestimates:
    """Storage class for the original parameter estimations
       CONTENTS:
         period [Days] - Initial orbital period
                         ***By default this is fixed during fitting***
         epoch [Days] - Initial epoch
         duration [Hours] - Initial duration fitted **In Hours**
         depth [ppm] - Initial depth
    """
    def __init__(self):
        self.period = 1.0
        self.epoch = 0.1
        self.duration = 3.0
        self.depth = 100.0

    def __str__(self):
        for k in self.__dict__:
            print k, self.__dict__[k]
        return ''

class trp_planetestimates:
    """Storage class for estimating a planet model based
       upon the trapezoid fit solution.  See Carter et al. 2008
       CONTENTS:
         u1 - quadratic limb darkening parameters to use
         u2 - ''
         period [Days] - Resulting period currently not fit
         radiusRatio - from purely geometric depth=(Rp/Rstar)^2
         impactParameter - Impact parameter
         tauzero - [Days] - transit timescale constant
         semiaxisRatio - Semi-major axis to stellar radius ratio
         surfaceBright - Limb darkened surface brightness at crossing
                         impact parameter
         equivRadiusRatio - Crude approximation to radius ratio
                            taking into account limb darkening
                            that works better than the purely geometric
                            radius ratio
         minDepth [ppm] - minimum depth from model
         avgDepth [ppm] - mean depth across transit
         epoch - epoch of fit midpoint
         bigT [day] - trapezoid model full duration
         littleT [day] - trapezoid model ingress/egress duration
         depth [ppm] - trapezoid model depth parameter
    """
    def __init__(self):
        self.u1 = 0.40 # limb darkening for Sun in Kepler passband
        self.u2 = 0.27
        self.period = 1.0
        self.radiusRatio = 0.0
        self.impactParameter = 0.5
        self.tauzero = 0.1
        self.semiaxisRatio = 20.0
        self.surfaceBright = 0.5
        self.equivRadiusRatio = 0.0
        self.minDepth = 0.0
        self.epoch = 0.0
        self.bigT = 0.0
        self.littleT = 0.0
        self.depth = 0.0
    def __str__(self):
        for k in self.__dict__:
            print k, self.__dict__[k]
        return ''

class trp_ioblk:
    """Define a class that contains all the data needed to perform
       trapezod fits.  Numerous functions will use as input this
       class.  This is purely a storage class
       See trp_setup to illustrate how to iinitialize these
       storage classes
       CONTENTS:
       parm - [class] Storage class trp_parameters for algorithm parameters
       origests - [class] Storage class trp_originalestimates for initial
                          parameter estimates
    """
    def __init__(self):
        self.parm = trp_parameters()
        self.origests = trp_originalestimates()
        self.planetests = trp_planetestimates()
        self.physval_names = ['']
        self.fixed = np.array([0])
        self.nparm = 0
        self.physval_mins = np.array([0.0])
        self.physval_maxs = np.array([0.0])
        self.physvals = np.array([0.0])
        self.physvalsavs = np.array([0.0])
        self.bestphysvals = np.array([0.0])
        self.boundedvals = np.array([0.0])
        self.boundedvalsavs = np.array([0.0])
        self.bestboundedvals = np.array([0.0])
        self.model = np.array([0.0])
        self.errscl = 1.0
        self.chi2min = 0.0
        self.minimized = False
        self.sampleit = np.array([0.0])
        self.fitdata = np.array(0, dtype=np.bool)
        self.normlc = np.array([0.0])
        self.normes = np.array([0.0])
        self.normts = np.array([0.0])
        self.normots = np.array([0.0])
        self.timezpt = 0.0

    def __str__(self):
        for k in self.__dict__:
            print k, self.__dict__[k]
        return ''

def boundedvals(ioblk):
    """Convert parameters to bounded versions that the minimzer will use
       INPUT:
         ioblk - [class] trp_ioblk class
       OUTPUT:
         ioblk - [class]
         err - [0 ok ; 1 not ok]
    """
    err = 0 # Error flag
    maxmindelta = ioblk.physval_maxs - ioblk.physval_mins
    datamindelta = ioblk.physvals - ioblk.physval_mins
    ioblk.boundedvals = -np.log( maxmindelta / datamindelta - 1.0)
    if ~np.isfinite(ioblk.boundedvals).all():
        print "Bounded Vals Bad"
        print ioblk.boundedvals
        print ioblk.physvals
        err = 1
    return ioblk, err

def unboundedvals(ioblk):
    """Convert bounded parameter values that the minimizer uses to physvals
       INPUT:
         ioblk - [class] trp_ioblk class
       OUTPUT:
         ioblk - [class]
         err - [0 ok ; 1 not ok]
    """
    err = 0 # Error flag
    maxmindelta = ioblk.physval_maxs - ioblk.physval_mins
    ioblk.physvals = ioblk.physval_mins + \
                     (maxmindelta / (1.0 + np.exp( -ioblk.boundedvals )))
    #if np.sum( np.isfinite(ioblk.physvals) ) != np.size(ioblk.boundedvals) :
    if ~np.isfinite(ioblk.physvals).all():
        print "UnBounded Vals Bad"
        print ioblk.boundedvals
        print ioblk.physvals
        err = 1
    return ioblk, err

def trapezoid(t, depth, bigT, littleT):
    """Trapezoid shape for model
       INPUT:
       t -  [float] vector of independent values to evaluate
                    trapezoid model
       depth - [float] depth of trapezoid
       bigT - [float] full trapezoid duration
       littleT - [float] 'ingress/egress' duration
       OUTPUT:
       output - [float] vector of trapezoid model values
    """
    output = np.full_like(t, 1.0)
    t = np.abs(t)
    output = np.where(t <= bigT/2.0 - littleT/2.0, 1.0 - depth, output)
    output = np.where(np.logical_and(t > bigT/2.0 - littleT/2.0, \
                      t < bigT/2.0 + littleT/2.0),  \
                      1.0 - depth + ((depth/littleT)* \
                      (t-bigT/2.0 + littleT/2.0)), output)
    return output

def trapezoid_model_onemodel(ts, period, epoch, depth, bigT, littleT, subsamplen):
    """Make a trapezoid model at the given input parameters.  This routine
        generates the ioblk class which is used in the transit model.
        You can save time if you want to generate many models by
        calling this function once to generate the ioblk and then call
        trapezoid_model_raw() to generate the models at other inputs
            bypassing some of the setup routines in this function.
       INPUT:
           ts - Mid cadence time stamps
           period - Period of signal ***assumed fixed during model generation**
           epoch - Estimated epoch of signal.  Must be on same system
                           as ts
           depth [ppm] - Model depth
           bigT [hr] -full transit duration in hours
           littleT [hr] - ingress time in hours
           subsample - Subsample each cadence by this factor
        OUTPUT:
            ioblk - structure class containing model ligh curve
                    located at ioblk.modellc
    """
    # Instantiate trp_ioblk class and fill in values
    ioblk = trp_ioblk()
    ioblk.parm.debugLevel = 0
    ioblk.parm.samplen = subsamplen
    ioblk.normots = ts
    ioblk.origests.period = period
    ioblk.origests.epoch = epoch
    ioblk.origests.depth = depth
    ioblk.origests.duration = bigT
    # Calculate this from timeSeries
    ioblk.parm.cadlen = np.median(np.diff(ts))
    ioblk = trp_setup(ioblk)
    # update the tratio
    ioblk.physvals[3] = littleT / bigT
    ioblk, err = boundedvals(ioblk)

    ioblk.physvalsavs = ioblk.physvals
    ioblk.boundedvalsavs = ioblk.boundedvals
    ioblk, err = trapezoid_model(ioblk)
    return ioblk

def trapezoid_model_raw(ioblk, epoch, depth, bigT, littleT):
    """If you have a preexisting ioblk from fit or trapezoid_model_onemodel()
        You can just call this function to get another model
        at a different epoch depth duration and ingress time
        ****period is not variable at this point call
        trapezoid_model_onemodel() instead
       INPUT:
           ioblk - pre-existing ioblk from fitting or trapezoid_model_onemodel()
           epoch - Estimated epoch of signal.  Must be on same system
                           as ts
           depth [ppm] - Model depth
           bigT [hr] -full transit duration in hours
           littleT [hr] - ingress time in hour
        OUTPUT:
            ioblk - structure class containing model ligh curve
                    located at ioblk.modellc
    """
    ioblk.physvals[0] = epoch - ioblk.origests.epoch
    ioblk.physvals[1] = depth / 1.0e6
    ioblk.physvals[2] = bigT / 24.0
    ioblk.physvals[3] = littleT / bigT
    ioblk, err = boundedvals(ioblk)

    ioblk, err = trapezoid_model(ioblk)
    return ioblk


def trapezoid_model(ioblk):
    """Generate a subsampled model at the current parameters
       INPUT:
       ioblk - [class] trp_ioblk class structure
       OUTPUT:
       ioblk - [class] modified ioblk
       err - [0 ok; 1 not ok] Error flag
    """
    err = 0
    to = ioblk.physvals[0]
    depth = ioblk.physvals[1]
    bigT = ioblk.physvals[2]
    tRatio = ioblk.physvals[3]
    littleT = tRatio * bigT
    per = ioblk.origests.period
    ts = ioblk.normts
    phi = phaseData(ts, per, to)
    lc = np.ones_like(ioblk.normts)
    cadlen = ioblk.parm.cadlen
    samplen = ioblk.parm.samplen
    # Call trapezoid model for data points without any subsampling needed
    idx = np.where(np.logical_and(ioblk.fitdata, ioblk.sampleit == 1))[0]
    if idx.size > 0:
        ztmp = phi[idx] * per
        lctmp = trapezoid(ztmp, depth, bigT, littleT)
        lc[idx] = lctmp
    # Call trapezoid model for data points that need subsampling
    idx = np.where(np.logical_and(ioblk.fitdata, ioblk.sampleit > 1))[0]
    if idx.size > 0:
        ztmp = phi[idx] * per
        deltaXSmall = cadlen / np.float(samplen)
        smallBlock = np.linspace(-cadlen/2.0 + deltaXSmall/2.0,
                                  cadlen/2.0 - deltaXSmall/2.0, samplen)
        oN = ztmp.size
        ztmp_highres = np.tile(ztmp, samplen)
        ztmp_highres = np.reshape(ztmp_highres, (samplen, oN))
        smallBlock_highres = np.tile(smallBlock, oN)
        smallBlock_highres = np.reshape(smallBlock_highres, (oN, samplen))
        smallBlock_highres = np.transpose(smallBlock_highres)
        ztmp_highres = ztmp_highres + smallBlock_highres
        ztmp_highres = ztmp_highres.ravel(order='F')
        lctmp_highres = trapezoid(ztmp_highres, depth, bigT, littleT)
        nN = ztmp_highres.size
        lctmp = lctmp_highres.reshape([oN, nN/oN]).mean(1)
        lc[idx] = lctmp

    ioblk.modellc = lc
    if np.sum(np.isfinite(lc)) != lc.size:
        err = 1
    return ioblk, err

def trp_setup(ioblk):
    """Setup various data products before minimizing
       INPUT:
       ioblk - [class] trp_ioblk class structure
       OUTPUT:
       ioblk - [class] modified ioblk
    """
    per = ioblk.origests.period
    eph = ioblk.origests.epoch
    dur = ioblk.origests.duration
    depth = ioblk.origests.depth / 1.0e6
    durday = dur / 24.0
    phidur = dur / 24.0 / per

    # Normalize the time series
    ts = ioblk.normots
    medianEvent = np.median(np.round((ts - eph)/per))
    ioblk.timezpt = eph + (medianEvent * per)
    ioblk.normts = ioblk.normots - ioblk.timezpt
    # identify in transit data to over sample and fitting region
    phi = phaseData(ioblk.normts, per, 0.0)
    ioblk.sampleit = np.where(abs(phi) < (phidur * 1.5), ioblk.parm.samplen, 1)
    ioblk.fitdata = np.where(abs(phi) < (phidur * ioblk.parm.fitregion),\
                             True, False)
    # always fit less than a 0.25 of phase space for stability
    #  and efficiency reasons
    ioblk.fitdata = np.where(abs(phi) > 0.25, False, ioblk.fitdata)

    # Set parameters and bounds
    ioblk.physval_names = ['To', 'Depth', 'BigT', 'TRatio']
    ioblk.physval_mins = np.array([-durday*1.5, 1.0e-6, 0.0, 1.0e-10])
    ioblk.physval_maxs = np.array([ durday*1.5, depth*5.0, durday*3.0, 1.0])
    ioblk.fixed = np.array([0, 0, 0, 0])
    ioblk.physvals = np.array([0.0, depth, durday, 0.2])
    ioblk.nparm = np.size(ioblk.fixed)

    ioblk.modellc = np.full_like(ioblk.normlc, 1.0)
    ioblk.chi2min = ioblk.normlc.size * 2000.0
    ioblk.likecount = 0
    ioblk, err = boundedvals(ioblk)

    # physvalsavs and boundedvalsavs are used to store parameters
    #  that are fixed during the calculation
    #  ***They must be populated with fixed values before moving forward
    ioblk.physvalsavs = ioblk.physvals
    ioblk.boundedvalsavs = ioblk.boundedvals

    ioblk.bestphysvals = ioblk.physvals
    ioblk.bestboundedvals = ioblk.boundedvals
    ioblk.minimized = False
    return ioblk

def trp_likehood(pars,ioblk):
    """Return a residual time series of data minus model
       trp_setup(ioblk) should be called before this function is called
       INPUT:
       pars - [numpy array] vector of parameter values
       ioblk - [class] trp_ioblk class structure
       OUTPUT:
       residuals - sum of squares of residuals of data - model
       ioblk - [class] modified ioblk
    """
    ioblk.likecount += 1
    # Update parameters into bounded values
    idx = np.where(ioblk.fixed == 0)[0]
    ioblk.boundedvals[idx] = pars
    ioblk.boundedvals = np.where(ioblk.fixed == 1, ioblk.boundedvalsavs,
                                 ioblk.boundedvals)
    # Convert to unbounded values
    ioblk, err = unboundedvals(ioblk)
    # Generate Model
    ioblk, err = trapezoid_model(ioblk)
    # Calculate residuals
    idx = np.where(ioblk.fitdata)[0]
    residuals = (ioblk.normlc[idx] - ioblk.modellc[idx])/(ioblk.normes[idx] * ioblk.errscl)
    # Return scalar summed residuals
    residuals = np.sum(residuals**2)

    # Do plotting
    if ioblk.parm.debugLevel > 2:
        if ioblk.likecount == 1: # Setup  figures for first time
            ioblk.fighandle = plt.figure(figsize=(3,2),dpi=300,
                                         facecolor='white')
            ioblk.axhandle = plt.gca()
            ioblk.axhandle.set_position([0.125, 0.125, 0.825, 0.825])
            ioblk.axhandle.set_axis_bgcolor('white')
        if np.mod(ioblk.likecount, ioblk.parm.likehoodmoddisplay) == 0 \
              or ioblk.likecount == 1:
            plt.figure(ioblk.fighandle.number)
            plt.cla()
            period = ioblk.origests.period
            tzero = ioblk.physvals[0]
            ts = ioblk.normts
            phi = phaseData(ts, period, tzero)
            plt.plot(phi,ioblk.normlc,'.',markersize=0.6)
            plt.plot(phi,ioblk.modellc,'.r',markersize=0.6)
            plt.pause(0.0001) # This line forces a draw it seems
                              # getting matplotlib to plot in a non blocking
                              # manner has a storied history on the web
                              # this method may fail in later versions
            if ioblk.parm.debugLevel > 3:
                raw_input("Press [ENTER]")
    return residuals

def trp_iterate_solution(ioblk, nIter):
    """Peform multiple iterations starting from random initial conditions
       return the best solution in a chi2 sense among the nIter iterations
    """
    bestChi2s = np.zeros(nIter)
    bestParameters = np.zeros((ioblk.physvals.size, nIter))
    gdFits = np.zeros(nIter, dtype=np.bool)
    depth = ioblk.origests.depth / 1.0e6
    for i in range(nIter):
        ioblk.physvals = ioblk.physval_mins + \
                         np.random.rand(ioblk.physvals.size) * \
                         (ioblk.physval_maxs - ioblk.physval_mins)
        # Force depth parameter to start at minimum half the depth
        if ioblk.physvals[1] < np.abs(depth/2.0):
            ioblk.physvals[1] = depth / 2.0
        # Replace random starts with parameters values that are fixed
        ioblk.physvals = np.where(ioblk.fixed == 1, ioblk.physvalsavs, \
                                    ioblk.physvals)
        ioblk, err = boundedvals(ioblk)
        freeidx = np.where(ioblk.fixed == 0)[0]
        startParameters = ioblk.boundedvals[freeidx]
        #usemethod = 'Nelder-Mead'
        usemethod = 'Powell'
        useoptions = {'xtol': 1e-5, 'ftol': 1e-5, 'maxiter': 2000, 'maxfev': 2000}
        #usemethod = 'CG'
        #useoptions = {'gtol': 1e-5, 'maxiter': 2000}
        allOutput = opt.minimize(trp_likehood, startParameters, args=(ioblk,), \
                                 method=usemethod, options=useoptions)
        ioblk.boundedvals[freeidx] = allOutput['x']
        ioblk.boundedvals = np.where(ioblk.fixed == 1, ioblk.boundedvalsavs, \
                                     ioblk.boundedvals)
        ioblk, err = unboundedvals(ioblk)
        chi2min = allOutput['fun']
        if ioblk.parm.debugLevel > 0:
            strout = "%s %d %s %f" % ("It: ",i," Chi2: ",chi2min)
            print strout
            print ioblk.physvals
        if np.isfinite(ioblk.physvals).all():
            gdFits[i] = True
            bestChi2s[i] = chi2min
            bestParameters[:,i] = ioblk.physvals

    # Done with iterations find the best one by chi2min
    bestMaskedIdx = np.argmin(bestChi2s[gdFits])
    ioblk.chi2min = bestChi2s[gdFits][bestMaskedIdx]
    ioblk.bestphysvals = bestParameters[:,gdFits][:,bestMaskedIdx]
    ioblk.physvals = ioblk.bestphysvals
    ioblk, err = boundedvals(ioblk)
    ioblk.bestboundedvals = ioblk.boundedvals
    if ioblk.parm.debugLevel > 0:
        strout = "%s %f" % ("Overall Best Chi2 Min: ",ioblk.chi2min)
        print strout
        print ioblk.physvals
    ioblk.minimized = True
    return ioblk

def trp_estimate_planet(ioblk):
    """Convert the trapezoid fit solution into a crude estimate
       of a planet model that is close to trapezoid solution
       This fills out values in trp_planetestimates class
    """
    if not ioblk.minimized:
        strout = "Warning getting planet estimates for non converged \
                  trapezoid fit.  Do not trust results"
        print strout
    ioblk.planetests.period = ioblk.origests.period
    ioblk.planetests.epoch = ioblk.timezpt + ioblk.bestphysvals[0]
    ioblk.planetests.bigT = ioblk.bestphysvals[2]
    ioblk.planetests.littleT = ioblk.bestphysvals[3] * \
                                    ioblk.planetests.bigT
    ioblk.planetests.depth = ioblk.bestphysvals[1]
    # call likehood to get best transit model
    idx = np.where(ioblk.fixed == 0)[0]
    resids = trp_likehood(ioblk.bestboundedvals[idx], ioblk)
    trapmodlc = ioblk.modellc
    ioblk.planetests.minDepth = (1.0 - trapmodlc.min()) * 1.0e6
    ioblk.planetests.radiusRatio = np.sqrt(ioblk.planetests.minDepth / 1.0e6)
    ioblk.planetests.impactParameter = np.sqrt(1.0 - \
                    np.amin([ioblk.planetests.radiusRatio * \
                    ioblk.planetests.bigT/ioblk.planetests.littleT, 1.0]))
    ioblk.planetests.tauzero = np.sqrt(ioblk.planetests.bigT * \
                    ioblk.planetests.littleT / 4.0 / \
                    ioblk.planetests.radiusRatio)
    ioblk.planetests.semiaxisRatio = ioblk.planetests.period / 2.0 / \
                    np.pi / ioblk.planetests.tauzero
    mu = np.sqrt(1.0 - ioblk.planetests.impactParameter**2)
    ioblk.planetests.surfaceBright = 1.0 - ioblk.planetests.u1*(1.0-mu) - \
                    ioblk.planetests.u2*(1.0-mu)**2
    ioblk.planetests.equivRadiusRatio = ioblk.planetests.radiusRatio / \
                    np.sqrt(ioblk.planetests.surfaceBright)

    return ioblk

def trapezoid_fit(timeSeries, dataSeries, errorSeries, \
                  signalPeriod, signalEpoch, signalDuration, signalDepth, \
                  fitTrialN=13, fitRegion=4.0, errorScale=1.0, debugLevel=0,
                  sampleN=15, showFitInterval=30):
    """Perform a trapezoid fit to a normalized flux time series
       Assumes all data has the same cadence duration
       Period is fixed during the trapezoid fit
       AUTHOR: Christopher J Burke
       INPUT:
           timeSeries - Mid cadence time stamps
           dataSeries - Normalized time series
           errorSeries - Error time series
           signalPeriod - Period of signal ***assumed fixed during model fit**
           signalEpoch - Estimated epoch of signal.  Must be on same system
                           as timeSeries
           signalDuration [hr] - Estimated signal duration ***In hours**
           signalDepth [ppm] - Estimated signal depth
           fitTrialN - How many trial fits to perform starting at random
                       initial locations.  Increase this if you find the
                       minimization is returning local minima
           fitRegion - Fit data within fitRegion*signalDuration of signalEpoch
           errorScale - Default 1.0 - Scale the errorbars by this factor
           debugLevel - 0 Show nothing; 1-Show some text about iterations
                        2 Show some more text; 3 - Show graphical fit in
                           progress; 4 - pause for each graphical fit
           sampleN - Subsample each cadence by this factor
           showFitInterval - If debugLevel >=3 the show every showFitInterval
                              function evaluation
        OUTPUT:
           ioblk - An instance of trp_ioblk which is a class used to store
                   all information pertaining to fit results
    """
    # Instantiate trp_ioblk class and fill in values
    ioblk = trp_ioblk()
    ioblk.parm.debugLevel = debugLevel
    ioblk.parm.samplen = sampleN
    ioblk.parm.likehoodmoddisplay = showFitInterval
    ioblk.fitregion = fitRegion
    ioblk.normlc = dataSeries
    ioblk.normes = errorSeries
    ioblk.errscl = errorScale
    ioblk.normots = timeSeries
    ioblk.origests.period = signalPeriod
    ioblk.origests.epoch = signalEpoch
    ioblk.origests.duration = signalDuration   # input duration is hours
    ioblk.origests.depth = signalDepth

    # Calculate this from timeSeries
    ioblk.parm.cadlen = np.median(np.diff(timeSeries))

    # setup some more variables
    ioblk = trp_setup(ioblk)
    # Find solution by trying random initial conditions
    ioblk = trp_iterate_solution(ioblk,fitTrialN)
    # Convert the trapezoid fit solution into a pseudo planet model parameters
    ioblk = trp_estimate_planet(ioblk)
    return ioblk

# Run the test of a trapezoid model fit in gaussian noise
if __name__ == "__main__":
    # Make some fake data
    dataSpan = 80.0 # in Days
    exposureLength = 1.0/48.0 # in Days simulating 48 cadences per day
    nData = dataSpan / exposureLength
    noiseLevel = 40.0 # noise per observation in ppm
    signalDepth = 300.0 # signal depth in ppm
    signalDuration = 5.0 / 24.0 # in Days
    signalDurationHours = signalDuration * 24.0
    signalPeriod = 10.4203 # in Days
    signalEpoch = 5.1 # in Days
    timeSeries = np.linspace(0.0, dataSpan, nData);
    dataSeries = 1.0 + np.random.randn(nData) / 1.0e6 * noiseLevel
    errorSeries = np.full_like(dataSeries,noiseLevel/1.0e6)
    # Instantiate trp_ioblk class and fill in values
    ioblk = trp_ioblk()
    ioblk.parm.samplen = 15
    ioblk.parm.cadlen = exposureLength
    ioblk.fitregion = 4.0
    ioblk.normlc = dataSeries
    ioblk.normes = errorSeries
    ioblk.normots = timeSeries
    ioblk.origests.period = signalPeriod
    ioblk.origests.epoch = signalEpoch
    ioblk.origests.duration = signalDurationHours # input duration is hours
    ioblk.origests.depth = signalDepth
    # setup some more variables
    ioblk = trp_setup(ioblk)
    ioblk.physvals = np.array([0.0, signalDepth/1.0e6, signalDuration, 0.1])
    # Make a model trapezoid light curve
    ioblk, err = trapezoid_model(ioblk)

    #Phase data
    phasedSeries = phaseData(timeSeries, signalPeriod, signalEpoch)
    # Insert signal
    phaseDuration = signalDuration / signalPeriod
    dataSeries = dataSeries * ioblk.modellc
    #plt.plot(phasedSeries, dataSeries, '.')
    #plt.show()
    #plt.plot(timeSeries, dataSeries, '.')
    #plt.show()

    # Test fitting
    ioblk = trapezoid_fit(timeSeries, dataSeries, errorSeries, \
                  signalPeriod, signalEpoch+0.001, signalDurationHours*0.9, \
                  signalDepth*1.1, \
                  fitTrialN=2, fitRegion=4.0, errorScale=1.0, debugLevel=3,
                  sampleN=15, showFitInterval=30)
    print ioblk
    # test generating model
    newioblk = trapezoid_model_onemodel(timeSeries, signalPeriod, \
                    signalEpoch, signalDepth, signalDurationHours, \
                    signalDurationHours*0.1, ioblk.parm.samplen)
    plt.close('all')
    plt.plot(phasedSeries, newioblk.modellc,'.b')
    newioblk = trapezoid_model_raw(newioblk, signalEpoch+0.05, signalDepth*1.5, \
                    signalDurationHours*2.0, signalDurationHours*2.0*0.2)
    plt.plot(phasedSeries, newioblk.modellc, '.r')
    plt.show()
