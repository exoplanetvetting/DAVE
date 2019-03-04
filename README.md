# DAVE™
**D**iscovery **A**nd **V**etting of K2 **E**xoplanets™

## Summary

This repository implements a pipeline to find and vet planets planets
using data from [NASA's K2](http://keplerscience.arc.nasa.gov), and [TESS](https://www.nasa.gov/tess-transiting-exoplanet-survey-satellite/) missions.

The pipeline performs the following steps:

1. Create a subset of targets for testing
2. Light Curve Generation, cotrend, detrend. For K2, removing instrumental effects e.g. PDC light curves Dan Foreman-Mackey's method; For TESS using Savitzky-Golay filter
3. Search for planets -- BLS
4. Fitting a planetary model
5. Produce vetting Metrics
6. Output useful data products


## Installation

*dave* depends on a few packages including numpy, scipy, matplotlib, and octave.
To check if you have all the dependencies, run:
```
python checkRequirements.py
```

You will almost certainly need to install @dfm's [python-bls](https://github.com/dfm/python-bls) package.


## Example use for K2
1. from dave.pipeline import main
2. cfg=main.loadMyConfiguration()
3. clip = main.runOne(206103150,cfg)

## Example use for TESS
1. from dave.tessPipeline import vet_tess_
2. detrendType_ = "eleanor"
3. clip = vet_tess_.runOneDv(detrendType, 2, 307210830, 1, 3.690613, 1356.2038, 1863, 1.27)
4. vet_tess_.runExport(clip,"tmp.txt")


where input for runOneDv is:
Sector, TIC ID, Planet Number, Period, BTJD, Transit Depth [ppm], Transit Duration [hours]. 

Currently supported "detrendType" are "tess" and "eleanor", where "tess" refers to the [SPOC short-cadence data](https://archive.stsci.edu/prepds/tess-data-alerts/), and "eleanor" refers to [eleanor-generated FFI lightcurves](http://adina.feinste.in/eleanor/)

## Additional requirements
To get DAVE to work for K2, you will have to create a directory called `.mastio/k2` in your home directory.

For K2, you will also need to download and extract
[kplr2011265_prf.tar.gz](https://archive.stsci.edu/pub/kepler/fpc/kplr2011265_prf.tar.gz)
into `.mastio/keplerprf`. 

For TESS, you will also need to download the TESS PRF Matlab files (https://archive.stsci.edu/tess/all_products.html). For more information, see designnote.pdf (Copyright Fergal Mullally).
