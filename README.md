# DAVE™
**D**iscovery **A**nd **V**etting of K2 **E**xoplanets™ ([Kostov, V. B., Mullally, S. E., Quintana , E. V. et al. 2019, AJ, 157, 3](https://arxiv.org/abs/1901.07459))

## Summary

This repository implements a pipeline to find and vet planets planets
using data from [NASA's K2](http://keplerscience.arc.nasa.gov), and [TESS](https://www.nasa.gov/tess-transiting-exoplanet-survey-satellite/) missions.

The pipeline performs the following steps:

1. Create a subset of targets for testing
2. Light Curve Generation, cotrend, detrend. For K2, removing instrumental effects using Savitzky-Golay filter
3. Search for planets -- BLS
4. Fitting a planetary model
5. Produce vetting Metrics
6. Output useful data products


## Installation

*dave* depends on a few packages including numpy, scipy, matplotlib, and octave. For K2, you will also need to install @dfm's [python-bls](https://github.com/dfm/python-bls) package; For TESS, you will need the [astropy BLS package](http://docs.astropy.org/en/stable/api/astropy.stats.BoxLeastSquares.html). To check if you have all the dependencies, run:
```
python checkRequirements.py
```
You may also need to pip install the following packages: astropy, metric_learn, sklearn ,pyfits, astroquery, conda install gnuplot (for which you will need Conda-Forge), PyGnuplot, parmap, clipboard, lpproj, and numba. If you are using Python 3+, and get an error when trying to run "from sklearn.utils.extmath import pinvh" in python, you will need to replace this line in file "~/site-packages/metric_learn/sdml.py" with "from scipy.linalg import pinvh".

## Example use for K2
from python:
1. from dave.pipeline import main
2. cfg=main.loadMyConfiguration()
3. clip = main.runOne(206103150,cfg)

## Example use for TESS
from python:
1. from dave.tessPipeline import vet_tess_
2. detrendType = "tess_2min"
3. clip = vet_tess_.runOneDv(detrendType, 2, 307210830, 1, 3.690613, 1356.2038, 1863, 1.27)
4. outfile = "tmp.txt"
5. export_ = vet_tess_.runExport(clip,outfile_)

where input for runOneDv is:
Sector, TIC ID, Planet Number, Period, BTJD, Transit Depth [ppm], Transit Duration [hours]. 

or from the terminal: "./runOneTESS"

Currently supported "detrendType" are "tess_2min" and "eleanor", where "tess_2min" refers to the [SPOC short-cadence data](https://archive.stsci.edu/prepds/tess-data-alerts/), and "eleanor" refers to [eleanor-generated FFI lightcurves](http://adina.feinste.in/eleanor/)

## Additional requirements
To get DAVE to work for K2, you will have to create a directory called `.mastio/k2` in your home directory.

For K2, you will also need to download and add
[Kepler's PRF files](https://archive.stsci.edu/missions/kepler/fpc/prf/)
into `.mastio/keplerprf`. 

For TESS, you will also need to download the TESS PRF Matlab files (https://archive.stsci.edu/tess/all_products.html). For more information, see designnote.pdf (Copyright Fergal Mullally).
