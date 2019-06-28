
import matplotlib.pyplot as mp
import numpy as np

from matplotlib.patches import Ellipse
from scipy.stats import chi2




def example():
    n = 100
    x = np.arange(n)
    y = x + .3*n*np.random.randn(n)


    mp.clf()
    mp.plot(x, y, 'bo')
    mp.plot(0,0, 'ko', ms=12)

    plotErrorEllipse(x, y)
    print(computeProbabilityOfObservedOffset(x, y))
    mp.grid()
    mp.axis([-15, 105, -60, 160])


def computeProbabilityOfObservedOffset(x, y, p=None):
    """Compute probability that p is consistent with mean of distribution.

    For a 2 dimensional distribution of points, given by ``x`` and ``y``,
    compute the probability that the mean of the distribution is consistent
    with the input point ``p``

    Inputs:
    -------------
    x, y
        (float) Input values. Require that ``len(x) == len(y)``

    Optional Inputs:
    -----------------
    p
        (array of length 2) Point to test. Default is [0,0]

    Returns:
    -----------
    probOffset
        (float) Probability that input point is consistent with mean
        of distribution.
    chiSquare
        (float) The chi squared of the point. For highly inconsistent points
        the computed probability of offset flatlines at zero. The chisquare
        value can then be used to estimate the relative consistencies of
        different points.


    Notes:
    ---------
    See ``plotErrorEllipse`` for a description of the algorithm
    """

    if p is None:
        p = [0, 0]
    p = np.array(p)
    assert(len(p) == 2)

    assert(len(x) == len(y))
    if len(x) < 2:
        raise ValueError("Need at least two points to compute probability of offset")

    mu = np.array([np.mean(x), np.mean(y)])
    cov = np.cov(x, y) / np.sqrt(len(x))

    eigenVals, eigenVecs = np.linalg.eigh(cov)
    v1 = eigenVecs[:, 0]
    v2 = eigenVecs[:, 1]

    pDash = (p-mu)
    offset_pix = np.array([np.dot(pDash, v1), np.dot(pDash, v2)])
    sigma = np.sqrt(eigenVals)
    offset_sigma = offset_pix / sigma

    s = np.sum(offset_sigma**2)
    probOffset = chi2.cdf(s, 2)

    return probOffset, s


def plotErrorEllipse(x, y, prob=[.68, .95, .997], **kwargs):
    """Compute and plot error ellipses around the mean of a 2d distribution.

    Given two arrays, ``x`` and ``y`` where the values of each are drawn
    from a Gaussian distribution, but the values of ``y`` are correlated with
    the values of ``x``, compute and plot both the mean of the distribution
    and the uncertainty ellipses around that mean.

    Inputs:
    -------------
    x, y
        (float) Input values. Require that ``len(x) == len(y)``

    Optional Inputs:
    -----------------
    prob
        (list or array) List of probability contours to draw. The defaults
        are 68%, 95% and 99.7%. Each value must be in the range (0,1)

    marker
        (string) Marker type to plot the location of the mean. Default: 'o'
    mfc
        (string) Colour of above marker. Default is 'r', can be any legal
        matplotlib colour specification.
    ms
        (int) Marker size. Default 12

    Any other optional inputs are passed to matplotlib.patches.Ellipse

    Returns:
    ----------
    **None**


    Notes:
    ---------
    Adapted from `www.visiondummy.com  <http://www.visiondummy.com/2014/04/draw-error-ellipse-representing-covariance-matrix/#Source_Code>`_

    The general idea is as follows:

    * Compute the covariance matrix for x and y.
    * Compute the eigen vectors of the covariance matrix, corresponding
      to vectors parallel to the major and minor axes.
    * The eigen values correspond to the variances of the distribution
      along the major and minor axes.
    * To compute the lengths of the major and minor axes for a given confidence
      interval, compute the chi^2 value that corresponds to that confidence
      interval. The length of the semi major axis is then
      :math:`\sqrt{ \chi^2 \sigma^2}`, where :math:`\sigma^2` is
      the variance.

    The above proceedure will give you confidence intervals for the
    distribution (e.g the p=.95 ellipse will contain 95% of the points in the
    distribution. To get the confidence interval on the mean, divde the
    covariance matrix in the first step by :math:`\sqrt{N}` where N is the
    number of input points.

    """

    assert(len(x) == len(y))
    prob = np.array(prob)
    assert(len(prob) >= 1)
    assert(np.all(prob > 0) & np.all(prob < 1))

    if len(x) < 2:
        raise ValueError("Need at least two points to generate centroid plots")

    # Default values for mp.plot
    marker = kwargs.pop('marker', 'o')
    mfc = kwargs.pop('mfc', 'r')
    ms = kwargs.pop('ms', 12)

    # Plot the mean position
    muX = np.mean(x)
    muY = np.mean(y)
    mp.plot(muX, muY, marker=marker, mfc=mfc, ms=ms)

    # Default values for ellipse
    if 'alpha' not in kwargs:
        kwargs['alpha'] = .2

    if 'color' not in kwargs:
        kwargs['color'] = 'r'

    # Dividing by root N gives the covariance matrix for the location of
    # the mean.
    cov = np.cov(x, y) / np.sqrt(len(x))

    eigenVals, eigenVecs = np.linalg.eig(cov)
    v1 = eigenVecs[:, 0]

    # Matplotlib's Ellipse takes an angle not a pair of vectors.
    angle_deg = np.degrees(np.arctan2(v1[1], v1[0]))

    ax = mp.gca()
    for p in prob:
        scale = inverseChiSquare(1-p)  # Convert prob to chisq
        sma = np.sqrt(scale * eigenVals[0])
        smi = np.sqrt(scale * eigenVals[1])
        ell = Ellipse(xy=[muX,muY], width=2 * sma, height=2 * smi,
                      angle=angle_deg, **kwargs)
        ax.add_patch(ell)


def inverseChiSquare(p, dof=2):
    """Compute inverse chi square function.

    For a given value of ``p``, computed the chi-square value
    corresponding to that probability.

    For example, the probability of getting a chi square value greater
    than 2.996 with two degrees of freedom is <0.05. So
    ``inverseChiSquare(.05, 2) == 2.996``

    Inputs:
    ------------
    p
        (float) Desired probability

    Optional Inputs
    ----------------
    dof
        (int) Degrees of freedom.


    Returns:
    -------------
    The chi-squared value corresponding to the input probability

    Notes:
    ---------
    I would expect scipy to support this, but it doesn't. The easiest
    way to compute it seems to be to interpolate over the availble
    chi-square distribution function.
    """

    y = np.logspace(2, -3, 100)
    chi = chi2(df=dof)
    x = 1 - chi.cdf(y)

    return np.interp(p, x, y)
