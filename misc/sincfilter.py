
import matplotlib.pyplot as mp
import numpy as np

"""
Functions to apply high and low pass filters to 1d data.
Based on "The Scientist and Engineers' Guide to Digital Signal
    Processing" by Steven Smith. Ch 16, Eqn 16.4. This book is available
    online  at www.dspguide.com/ch16/

Remember, these functions apply to evenly spaced data only.
"""

__version__ = "$Id: sincfilter.py 1983 2015-03-08 04:50:13Z fmullall $"
__URL__ = "$URL: svn+ssh://flux/home/fmullall/svn/kepler/k2phot/sincfilter.py $"


def normalisedFrequencyFromPeriod(period, timespan):
    return normaliseFrequency(1/float(period), timespan)

def normaliseFrequency(f, timespan):
    """Convert a frequency in the units of the data into normalised frequency.

    Useful for computing input arguments in {high|low}Pass()

    Inputs:
    f           (float) Frequency in, e.g, Hz
    timespace   (float) Range of data in , eg, seconds. This is the time
                interval between the last and first data points in your
                set.
    """
    return f/float(timespan)





def blackman(sysSize):
    """Generate the Blackman apodisation.

    Taken from "The Scientist and Engineers' Guide to Digital Signal
    Processing" by Steven Smith. Ch 16, Eqn 16.4. This book is available
    online  at www.dspguide.com/ch16/2.html

    This function is used to apodise a sinc filter to improve it's
    stopband attentuation.

    sysSize:    (int) Number of points being used in the filter. Smith
                calls this M

    Returns:
    A 1d array, blackman[i]. The sinc filter is apodised by computing
    apodisedFilter[i] = sinc[i] * blackman[i]
    """

    sysSize = float(2*sysSize+1)
    i = np.arange(sysSize)

    #Compute function from right to left. No reason other than it's
    #easier to implement
    arg = 4.*np.pi*i/sysSize
    blackman = .08*np.cos(arg)

    arg = 2.*np.pi*i/sysSize
    blackman -= 0.5*np.cos(arg)

    blackman += 0.42
    return blackman


def rectangle(sysSize):
    return np.ones(2*sysSize+1)





def highPass(y, normalisedCutoffFreq, filterHalfWidth, apodise=blackman):
    """Apply a high pass filter to evenly spaced data

    Based on sincFilter in this module. See that function for more
    detailed info on who Smith is.

    Inputs:
    y   (1d float array)    Array of evenly spaced data points to filter
    normalisedCutOffFreq    (float) Frequency at which 50% attentuation
                            occurs. This input must have a value between
                            [0,.5). Smith names this value f_c
    numPointsInFilter       (int)   Half the number of points included
                            in the filter. Smith refers to this value as M.
                            The number of points used is (2M+1)

    apodise                 (func) What apodisation function to use. The
                            default is usually the best choice. If you
                            don't want apodisation for some reason, set
                            this argument to rectangle

    Returns:
    A 1d array high pass filtered.
    """

    lp = lowPass(y, normalisedCutoffFreq, filterHalfWidth, apodise)
    return y - lp



def lowPass(y, normalisedCutoffFreq, filterHalfWidth, apodise=blackman):
    """Apply a low pass filter to evenly spaced data

    Based on sincFilter in this module. See that function for more
    detailed info on who Smith is.

    Inputs:
    y   (1d float array)    Array of evenly spaced data points to filter
    normalisedCutOffFreq    (float) Frequency at which 50% attentuation
                            occurs. This input must have a value between
                            [0,.5). Smith names this value f_c
    numPointsInFilter       (int)   Half the number of points included
                            in the filter. Smith refers to this value as M.
                            The number of points used is (2M+1)

    apodise                 (func) What apodisation function to use. The
                            default is usually the best choice. If you
                            don't want apodisation for some reason, set
                            this argument to rectangle

    Returns:
    A 1d array low pass filtered.
    """

    assert(filterHalfWidth %2 == 0)
    fc = normalisedCutoffFreq
    num = filterHalfWidth


    filt = apodise(num) * sincFilter(fc, num)
    filt /= np.sum(filt)

    #mp.subplot(211)
    #mp.plot(filt)

    #mp.subplot(212)
    #ft = np.fft.rfft(filt)/ np.sqrt(len(filt))
    #mp.semilogy(np.abs(ft))

    out = np.empty(len(y))
    for i in range(len(out)):
        i1 = max(i-num, 0)
        i2 = min(i+num, len(out)) - 1

        j1 = max(num-i, 0)
        j2 = min(j1 + i2 - i1, 2*num+1)

        out[i] = np.sum( y[i1:i2] * filt[j1:j2])


    return out





def sincFilter(normalisedCutoffFreq, numPointsInFilter):
    """
    Compute a sincFilter with numPointsInFilter points

    Taken from "The Scientist and Engineers' Guide to Digital Signal
    Processing" by Steven Smith. Ch 16, Eqn 16.4. This book is available
    online  at www.dspguide.com/ch16/2.html

    A sinc filter is a low pass filter with a sharp cutoff. The larger
    the value of numPointsInFilter, the sharper the cutoff.

    Sinc functions suffer from ringing on either side of the cutoff which
    can be dramatically reduced by apodising the filter with, e.g., a
    blackman filter.


    The
    Inputs:
    normalisedCutOffFreq    (float) Frequency at which 50% attentuation
                            occurs. This input must have a value between
                            [0,.5). Smith names this value f_c
    numPointsInFilter       (int)   Half the number of points included
                            in the filter. Smith refers to this value as M.
                            The number of points used is (2M+1)

    Returns:
    A 1d array of length numPointsInFilter.

    Notes:
    If your have a desired cutoff frequency in Hz, the normalised value
    can be computed as f/T, where f is the cutoff in Hz, and T is the
    timespan of the data in seconds.
    """

    if numPointsInFilter %2 != 0:
        raise ValueError("numPointsInFilter must be even")

    #Mnuemonics
    fc = normalisedCutoffFreq
    num = 2*numPointsInFilter + 1
    assert(fc >= 0)
    assert(fc < 0.5)

    i = np.arange(num)
    j = i- .5*num

    #return np.sinc( 2*np.pi*fc * j)
    numerator = np.sin(2*np.pi*fc * j)
    denom = np.pi*j

    #Catch the division by zero, and replace with 1
    idx = j == 0
    numerator[idx] = 1
    denom[idx] = 1
    return numerator/denom






def example():
    n =2e4
    fc = .05
    sysSize = 100

    mp.clf()

    y = np.random.randn(n)
    y1 = lowPass(y, fc, sysSize, apodise=rectangle)
    y2 = lowPass(y, fc, sysSize, apodise=blackman)

    ft = np.fft.rfft(y) / np.sqrt(n)
    #ft1 = np.fft.rfft(y1) / np.sqrt(n)
    #ft2 = np.fft.rfft(y2) / np.sqrt(n)

    mp.clf()
    #mp.subplot(211)
    #mp.plot(y, 'b.')
    ##mp.plot(y1, 'r.')
    #mp.plot(y2, 'g.')


    #mp.subplot(212)
    #mp.plot(np.abs(ft), 'b-')
    ##mp.plot(np.abs(ft1), 'r-')
    #mp.plot(np.abs(ft2), 'g-')



    fc = .01
    sysSize = 400
    y =  sincFilter(fc, sysSize)
    mp.plot(y, 'b-')
    y = blackman(sysSize) * sincFilter(fc, sysSize)
    mp.plot(y, 'g-')
