# dave™
**D**iscovery **A**nd **V**etting of K2 **E**xoplanets™

## Summary

This repository implements a pipeline to find and vet planets planets
using data from [NASA's K2 and TESS missions](http://keplerscience.arc.nasa.gov).

The pipeline performs the following steps:

1. Create a subset of targets for testing
2. Light Curve Genertion/CoTrend -- remove instrumental effects
	e.g. PDC light curves
	     Dan Foreman-Mackey's method
3. Detrend/Search for planets -- BLS
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
from dave.pipeline import main
cfg=main.loadMyConfiguration()
clip = main.runOne(206103150,cfg)

## Example use for TESS
Run "runOneTESS"

To get this to work, you will have to create a directory called `.mastio/k2` in your home directory.

For K2, you will also need to download and extract
[kplr2011265_prf.tar.gz](https://archive.stsci.edu/pub/kepler/fpc/kplr2011265_prf.tar.gz)
into `.mastio/keplerprf`. 

For TESS, you will need to download the TESS PRF Matlab files (https://archive.stsci.edu/tess/all_products.html). For more information, see designnote.pdf (Copyright Fergal Mullally).
