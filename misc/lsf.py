
import numpy as np

__version__ = "$Id: lsf.py 2209 2015-12-03 00:41:53Z fmullall $"
__URL__ = "$URL: svn+ssh://flux/home/fmullall/svn/kepler/py/lsf.py $"


##
#Common functions to fit
##
def poly(x, order, **kwargs):
    """Polynomial of the form \Sigma a_i x**i

       This is a working prototype of the analytic model to
       be used by lsf. This function is never called by the user
       directly.

       In general, the Lsf class fits an analytic function of
       the form f(x_i) = \sigma a_j g_j(x_i),
       where a_i is a fittable parameter, and g_j is a
       function of x_i. Calling poly with a given order
       returns the value of the derivative of f(x) with respect
       to a_order = g_order(x_i)

       g_j can depend on other parameters. These are passed
       in with kwargs.

       Inputs:
       x        (1d numpy array) ordinates
       order    (int) Which parameter to take the derivative with
                respect to
       kwargs   (dict) Optional extra arguments to pass to the
                function. poly doesn't require any.
    """


    return x**order


def sine(x, order, **kwargs):
    """Single sine curve of known period.

    Fits the equivalent of A.sin( 2 pi x /period)
    The equation is re-written as r1 sin(wx) + r2 cos(wx)
    where w = 2pi/period, and r1 and r2 are the fit coefficients

    To convert these parameters to the more useful amplitude
    and phase, see computeAmplitudeAndPhase below

    Period to fit should be passed in as kwarg with key 'period'

    """

    period = kwargs['period']
    w = 2*np.pi/period
    if order == 0:
        return np.sin(w * x)
    elif order == 1:
        return np.cos(w * x)
    else:
        raise ValueError("Order should be zero or one")


def fixedPowerLaw(x, order, **kwargs):
    """Fix a series of power laws where the exponent is known

    Example:
    A1 x**(b1) + A2 x**(b2) + ...
    where b_i is known, and A_i is to be fit.


    Inputs:
    x        (1d numpy array) ordinates
    order    (int) term of function to fit

    Required kwargs:
    exponents   (float or iterable) Exponents of power laws to fit (the
                b_i's above


    Notes:
    the optional argument exponents must be passed to Lsf when calling
    this function.

    order must equal number of exponents

    To fit a function where the exponents are not known you need a
    non-linear fit
    """


    exponents = kwargs.pop('exponents', None)
    if exponents is None:
        raise ValueError("fixedPowerLaw requires list of expoenents passed as kwarg")

    if not hasattr(exponents, "__len__"):
        exponents = [exponents]

    assert(order <= len(exponents))
    return x**exponents[order]



class Lsf():
    """Least squares fit to an analytic function

    Based on lsf.c in Wqed, which in turn is based on
    Bevington and Robinson.

    Numpy and scipy do least squares, but they don't make it easy
    to set weights, or get uncertainties in fit parameters
    """

    def __init__(self, x, y, s, order, func=poly, **kwargs):
        """
        Input
        x   (1d numpy array) ordinate (eg time)
        y   (1d numpy array) coordinate (eg flux)
        y   1 sigma uncertainties. Can be None, a scalar or
            a 1d numpy array

        order   (int) How many terms of func to fit
        func    (function object) The analytic function to fit
                see notes on poly() below. Default is poly
        kwargs  Additional arguments to pass to func.
        """

        if len(x) == 0:
            raise ValueError("x is zero length")

        if len(x) != len(y):
            raise IndexError("x and y must have same length")
        self.x = x
        self.y = y

        if order < 1:
            raise ValueError("Order must be at least 1")

        if order >= len(x):
            raise ValueError("Length of input must be at least one greater than order")
        size = len(x)
        if s is None:
            self.s = np.ones(len(x))
        elif np.isscalar(s):
            self.s = np.ones((size)) * s
        else:
            if len(x) != len(s):
                raise IndexError("x and s must have same length")
            self.s = s

        if not np.all(np.isfinite(self.x)):
            raise ValueError("Nan or Inf found in x")

        if not np.all(np.isfinite(self.y)):
            raise ValueError("Nan or Inf found in y")

        if not np.all(np.isfinite(self.s)):
            raise ValueError("Nan or Inf found in s")

        self.order = order
        self.func = func
        self.kwargs = kwargs

        self.fit()

    #def fit(self, x, y, s, order, func, **kwargs):
    def fit(self):
        """Fit the function to the data and return the best fit
        parameters"""
        dataSize = len(self.x)

        #Check array lengths agree

        df = np.empty((dataSize, self.order))
        for i in range(self.order):
            df[:,i] = self.func(self.x, i, **self.kwargs)
            df[:,i] /= self.s

        A = np.matrix(df.transpose()) * np.matrix(df)

        try:
            covar = A.I #Inverts
        except np.linalg.LinAlgError as e:
            print("lsf.py:Lsf.fit(): Error in fit: %s" %(e))
            self.param = None
            self.cov = None
            import pdb; pdb.set_trace()

        wy = np.matrix(self.y / self.s)
        beta = wy * df
        params = beta * covar


        #Store results and return
        self.param = np.array(params)[0]   #Convert from matrix
        self.cov = covar

        return self.param


    def getParams(self):
        return self.param


    def getVariance(self):
        """Variance is the sum of the squares of the residuals.
        If your points are gaussian distributed, the standard deviation
        of the gaussian is sqrt(var).
        """

        resid = self.getResiduals()

        val = np.sum( resid**2)
        val /= len(self.x) - self.order
        return val


    def getCovariance(self):
        """Get the covariance matrix

           Return:
           np.matrix object
        """
        assert(self.cov is not None)
        return self.cov

    def getUncertainty(self):
        """Get the uncertainty on the best fit params
            Return:
            1d numpy array
        """
        assert(self.cov is not None)
        unc = np.sqrt(self.cov.diagonal())
        return np.array(unc)[0]

    def getUncertaintyFromVariance(self):
        """Use variance of residuals to compute uncertainty

        Sometimes you have data with no easily accessible
        uncertainties. However, if you can assume
        a) The model is appropriate, and your residuals are pure noise
        b) each point has approx the same uncertainty

        you can use this function to compute the uncertainty
        by using the scatter in the residuals to estimate the
        average per data-point uncertainty.

        Inputs:
        (none)

        Returns:
        Vector of uncertainties
        """

        var = self.getVariance()
        return np.sqrt(var) * self.getUncertainty()


    def getChiSq(self):
        """Get chi square for fit

            Return double
        """
        chi = self.getResiduals() / self.s
        chisq = chi**2
        return np.sum(chisq)


    def getReducedChiSq(self):
        """Get reduced chi-squared."""
        num = len(self.x) - self.order
        return self.getChiSq()/float(num)


    def getBestFitModel(self, x=None):
        """Get best fit model to data

        Inputs:
        x   (1d numpy array) Ordinates on which to compute
            best fit model. If not specified use ordinates
            used in fit

        Return
        1d numpy array
        """

        assert(self.param is not None)

        if x is None:
            x = self.x

        try:
            len(x)
        except TypeError:
            x = np.asarray([x])

        par = self.param
        bestFit = np.zeros( (len(x)))
        for i in range(self.order):
            bestFit += par[i] * self.func(x, i, **self.kwargs)
        return bestFit

    def getResiduals(self):
        """Get residuals, defined as y - bestFit

        Returns:
        1d numpy array
        """
        return self.y - self.getBestFitModel()


    def getModelUncertainty(self, x=None):
        """Get the uncertainty in the value of the best fit model
        by propegating the uncertainties in the parameters and their
        covariances.

        For best results, use on data whose ordinates are centred on
        zero, or use a fitting function that does that translation for you.

        Inputs:
        ---------
        x:   Float or 1d array. Values to compute uncertainty at. Defaults
             to ordinates used in fit


        Returns:
        -----------
        A numpy array of the same length as x


        Caution:
        ----------
        Use on data where the ordinate is not centred on zero is not
        advised.

        This method is not well tested. Use it skeptically. It should
        work for polynomials, but I don't know how it will behave for
        more complicated functions.'

        Example:
        -------------
        >>> fObj = Lsf(x,y,s, order=2)
        >>> f = fObj.getBestModelFit()
        >>> df = fObj.getModelUncertainty()
        >>> mp.plot(x,y, 'ko')
        >>> mp.plot(x, f, 'r-')
        >>> mp.fill_between(x, f-df, f+df, 'k', alpha=.2)
        """
        assert(self.param is not None)

        if x is None:
            x = self.x

        try:
            len(x)
        except TypeError:
            x = np.asarray([x])


        unc = 0
        nPar = self.order
        for i in range(nPar):
            #cov[i,i] = sigma_i**2
            unc += self.cov[i,i] * self.func(x, i)**2
            for j in range(i+1, nPar):
                #cov[i,j] = sigma_ij not squarted
                tmp = self.cov[i,j] * self.func(x, i) * self.func(x, j)
                unc += 2*tmp**2

        return np.sqrt(unc)


def computeAmplitudeAndPhase(par):
    """Converts the results from fitting sine function to
    amplitude and phase"""

    amp = np.hypot(par[0], par[1])

    #r1 is cos(phase), r2 is -sin(phase) (because fitting wx-phi)
    # You can remove the - r2 if you add a minus to the lineAndSine()
    # for case 3. Don't do this, you'll only confuse yourself
    r1 = par[0]/amp
    r2 = -par[1]/amp

    #Remember the sign of the cos and sin components
    invcos = np.arccos(np.fabs(r1))
    invsin = np.arcsin(np.fabs(r2))

    #Decide what quadrant we're in
    #if(r1>=0 && r2 >=0) //1st quadrant, do nothing
    if r1 <=0 and r2>=0:
        #Second quadrant
        invcos = np.pi - invcos
        invsin = np.pi - invsin
    elif r1 <=0 and r2<=0:
        #Third quadrant
        invcos += np.pi
        invsin += np.pi
    elif r1>=0 and r2 <= 0:
        #4th quadrant
        invcos = 2*np.pi - invcos
        invsin = 2*np.pi - invsin


    #Choose the average of the two deteminations to
    #reduce effects of roundoff error
    phase = .5*(invcos+invsin)
    return amp, phase



def computeAmplitudeAndPhaseWithUnc(fObj):
    """Compute amplitude and phase errors for sine fit

    Taken from Appendix 1 of Breger (1999, A&A 349, 225)
    (said Appendix written by M Montgomery)

    Inputs
    ---------
    fObj
        Results of an least squares fit by the Lsf object of the `sine` function
        
    Returns
    -------
    A list:
    amplitude, phase, amplitude_unc, phase_unc
    """
    assert fObj.func == sine, "computeAmplitudeAndPhaseWithUnc only works for sine function"

    nPts = len(fObj.x)
    var = fObj.getVariance()
    amp, phase = computeAmplitudeAndPhase( fObj.getParams())
    
    ampUnc = np.sqrt( 2 * var / float(nPts) )     #A9 where \sigma(m) = sqrt(var)
    phaseUnc = ampUnc/amp                       #A14 (with some re-arranging)
    return amp, phase, ampUnc, phaseUnc
