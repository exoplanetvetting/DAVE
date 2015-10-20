from __future__ import division, print_function, absolute_import
import numpy as np

from .clean_and_search import Search

# call for bls search:
# outs = doSearch(time, flux, minPeriod, maxPeriod, )
#
# outputs:
#
# outs =
# {convolved_bls': array([ 0.00033225,  0.00034177,  0.00035353, ...,  0.0001447 ,
#          0.00013177,  0.00012505]),
#  'depth': 0.011876680987410373,
#  'duration': 5.4358004743484596,
#  'epoch': 2145.8181124895377,
#  'flux': array([ 0.00050952,  0.00061208,  0.0007578 , ..., -0.00015678,
#          0.        ,  0.00011196]),
#  'maxPeriod': 23.04694614343801,
#  'minPeriod': 0.5,
#  'period': 4.1591708630839328,
#  'bls_search_periods': array([ 23.04694614,  22.98218199,  22.9177808 , ...,   0.50009172,
#           0.50006114,   0.50003057]),
#  'time': array([ 2144.12322554,  2144.1436576 ,  2144.16408986, ...,  2213.2232006 ,
#          2213.24363228,  2213.26406397]),



def doSearch(time, flux, minPeriod, maxPeriod, ):

    period, epoch, duration, depth, bls_search_periods, \
            convolved_bls = do_ktwo_bls(time, flux, minPeriod, maxPeriod)

    return period, epoch, duration, depth, bls_search_periods, convolved_bls



def do_ktwo_bls(time, flux, min_period, max_period):
	S = Search(time, flux + 1, np.ones_like(flux)*1.E-5)
	S.do_bls2(min_period = min_period,
			  max_period = max_period,
			  min_duration_hours = 1.5,
			  max_duration_hours = 6.,
			  freq_step = 1.E-4,
			  doplot = False,
			  norm = False)


	# output from S:
	#    self.ingresses = ingresses
    #    self.egresses = egresses
    #    self.t_number = t_number
    #    self.epoch = epoch
    #    self.periods = periods
    #    self.bper = bper
    #    self.bpow = bpow
    #    self.depth = depth
    #    self.qtran = qtran
    #    self.duration = duration
    #    self.f_1 = f_1
    #    self.convolved_bls = convolved_bls
    #    self.approx_duration = approx_duration

	return S.periods[0], S.epoch, S.duration, S.depth, S.f_1, S.convolved_bls












