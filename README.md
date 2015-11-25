# dave™
**D**iscovery **A**nd **V**etting of K2 **E**xoplanets™

## Summary

This repository implements a pipeline to find and vet planets planets
using data from [NASA's K2 mission](http://keplerscience.arc.nasa.gov).

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


## Example use

```
from dave.pipeline import pipeline
cfg = pipeline.loadDefaultConfig()
pipeline.runOne(206103150, cfg)
```

(this will crash at present)

```
