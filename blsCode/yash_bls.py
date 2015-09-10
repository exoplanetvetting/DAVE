from __future__ import division, print_function
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import clean_and_search
from ktransit import FitTransit
from multiprocessing import Pool
from scipy import ndimage
import glob, timeit, sys
import time as pythonTime
# OPTIONS
doPlot = True
plotOption = 'save'
secondary = True
resultsFilename = '/Users/Yash/Desktop/results.txt'
figureSaveLocation = '/Users/Yash/Desktop/'
# -------- PLOTTING OPTIONS -------- #
import matplotlib

def plateau(array, threshold):
	"""Find plateaus in an array, i.e continuous regions that exceed threshold
	Given an array of numbers, return a 2d array such that
	out[:,0] marks the indices where the array crosses threshold from
	below, and out[:,1] marks the next time the array crosses that
	same threshold from below.
	Inputs:
	array       (1d numpy array)
	threshold   (float or array) If threshold is a single number, any point
				above that value is above threshold. If it's an array,
				it must have the same length as the first argument, and
				an array[i] > threshold[i] to be included as a plateau
	Returns:
	Numpy 2d array with 2 columns.
	Notes:
	To find the length of the plateaus, use
	out[:,1] - out[:,0]
	To find the length of the largest plateau, use
	np.max(out[:,1] - out[:,0])
	The algorithm fails if a value is exactly equal to the threshold.
	To guard against this, we add a very small amount to threshold
	to ensure floating point arithmetic prevents two numbers being
	exactly equal."""
	arr  = array.astype(np.float32)
	arr = arr - threshold + 1e-12
	arrPlus = np.roll(arr, 1)
	#Location of changes from -ve to +ve (or vice versa)
	#Last point is bogus , so we calculate it by hand
	sgnChange = arr*arrPlus
	#Roll around can't compute sign change for zeroth elt.
	sgnChange[0] = +1
	if arr[0] > 0:
		sgnChange[0] = -1
	loc = np.where(sgnChange < 0)[0]
	if np.fmod( len(loc), 2) != 0:
		loc.resize( (len(loc)+1))
		loc[-1] = len(arr)
	return loc


def outlierRemoval(time, flux):
	fluxDetrended = medianDetrend(flux, 3)
	out1 = plateau(fluxDetrended, 5 * np.std(fluxDetrended))
	out2 = plateau(-fluxDetrended, 5 * np.std(fluxDetrended))
	if out1 == [] and out2 == []:
		singleOutlierIndices = []
	else:
		outliers = np.append(out1, out2).reshape(-1,2)
		# Only want groups of one outlier, since > 1 may be transit points
		singleOutlierIndices = np.sort(outliers[(outliers[:,1] - outliers[:,0] == 1)][:,0])
	# Check periodicity of outliers, with PRECISION of 0.0205 days
	# 0.0205 days = 29.52 minutes = ~length of long cadence
	precision = 0.0205
	outlierTimes = time[singleOutlierIndices]
	diffs = [outlierTimes[i+1] - outlierTimes[i] for i in range(0, len(outlierTimes)-1)]
	diffs = [round(d, 5) for d in diffs]
	if len(singleOutlierIndices) >= 4:
		if len(set(diffs)) == len(diffs):
			possibleTimes = np.array([])
		else:
			period = max(set(diffs), key = diffs.count) # period = most common difference
			epoch = outlierTimes[diffs.index(period)]
			possibleTimes = np.arange(epoch, outlierTimes[-1] + 0.5*period, period)
		notOutliers = []
		for i in range(len(outlierTimes)):
			if np.any((abs(possibleTimes - outlierTimes[i]) < precision)):
				notOutliers.append(i)
		singleOutlierIndices = np.delete(singleOutlierIndices, notOutliers)
	elif len(singleOutlierIndices) == 3:
		if abs(diffs[0] - diffs[1]) < precision:
			singleOutlierIndices = []
	# Uncomment to see how the plotting algorithm worked for a lightcurve
	# ----------------------------- PLOTTING ----------------------------- #
	# plt.subplot(311)
	# plt.scatter(time, flux, marker = '.', s = 1, color = 'k', alpha = 1)
	# plt.scatter(time[singleOutlierIndices], flux[singleOutlierIndices], 
	#             s = 30, marker = 'o', facecolors = 'none', edgecolors = 'r')
	# plt.title('Original')
	# plt.subplot(312)
	# plt.scatter(time, fluxDetrended, marker = '.', s = 1, color = 'k', alpha = 1)
	# plt.scatter(time[singleOutlierIndices], fluxDetrended[singleOutlierIndices], 
	#             s = 30, marker = 'o', facecolors = 'none', edgecolors = 'r')
	# x1, x2, y1, y2 = plt.axis()
	# plt.hlines([-5*np.std(fluxDetrended), 5*np.std(fluxDetrended)], x1, x2, 
	#            color = 'b', linestyles = 'dashed')
	# plt.axis([x1, x2, y1, y2])
	# plt.title('Detrended')
	# plt.subplot(313)
	# plt.scatter(np.delete(time, singleOutlierIndices), np.delete(flux, singleOutlierIndices), 
	#             marker = '.', s = 1, color = 'k', alpha = 1)
	# plt.title('Outliers removed: ' + str(len(singleOutlierIndices)))
	# plt.show()
	# -------------------------------------------------------------------- #
	return np.delete(time, singleOutlierIndices), np.delete(flux, singleOutlierIndices)

def medianDetrend(flux, binWidth):
	halfNumPoints = binWidth // 2
	medians = []
	for i in range(len(flux)):
		if i < halfNumPoints:
			medians.append(np.median(flux[:i+halfNumPoints+1]))
		elif i > len(flux) - halfNumPoints - 1:
			medians.append(np.median(flux[i-halfNumPoints:]))
		else:
			medians.append(np.median(flux[i-halfNumPoints : i+halfNumPoints+1]))
	return flux - medians

def getPhase(time, flux, period, epoch, centerPhase = 0):
	"""Get the phase of a lightcurve.
	How it works using an example where epoch = 2, period = 3:
	1.	Subtract the epoch from all times [1, 2, 3, 4, 5, 6, 7] to get 
		[-1, 0, 1, 2, 3, 4, 5] then divide by the period [3] to get all time 
		values in phase values which gets you [-0.3, 0, 0.3, 0.6, 1, 1.3, 1.6]
	2.	Subtract the PHASE NUMBER (floor function) from each PHASE (date1)
		which gets you [0.7, 0, 0.3, 0.6, 0, 0.3, 0.6]
	3.	Sort all the adjusted phases to get [0, 0, 0.3, 0.3, 0.6, 0.6, 0.7]
		THERE WILL BE negative values in the beginning here, just not in this example
		since no ex. time value divided by the period left a decimal less than 0.25
	4.	Sort the flux values in the same way the phases were sorted
	Inputs:
	time			Time values of data. (IN DAYS)
	flux			Flux values of data.
	period 			Period of transit.
	epoch 			Epoch of transit.
	centerPhase 	Which phase should be at the center.
	Returns:
	q1				Phase values. (IN HOURS)
	f1				Flux values for each phase.
	"""
	epoch += centerPhase * period
	date1 = (time - epoch) / period  + 0.5
	phi1 = ((date1) - np.floor(date1)) - 0.5
	q1 = np.sort(phi1) * period * 24.
	f1 = flux[np.argsort(phi1)]
	return q1, f1

def fitModel(time, flux, guessDict, freeParPlanet, ferr = 0):
	if not np.all(ferr): ferr = np.ones_like(flux)*1.E-5
	freeParStar = ['rho']
	# Make the fitting object according to guess dictionary
	fitT = FitTransit()
	fitT.add_guess_star(ld1 = 0, ld2 = 0)
	fitT.add_guess_planet(period = guessDict['period'],
							  T0 = guessDict['T0'])
	fitT.add_data(time = time, flux = flux, ferr = ferr)
	fitT.free_parameters(freeParStar, freeParPlanet)
	fitT.do_fit()
	return fitT

def do_bls_and_fit(time, flux, min_period, max_period):
	S = clean_and_search.Search(time, flux + 1, np.ones_like(flux)*1.E-5)
	S.do_bls2(min_period = min_period,
			  max_period = max_period,
			  min_duration_hours = 1.5,
			  max_duration_hours = 6.,
			  freq_step = 1.E-4,
			  doplot = False,
			  norm = False)
	guessDict = {'period':   S.periods[0],
				 'T0':       S.epoch}
	freeParPlanet = ['period', 'T0', 'rprs']
	fitT = fitModel(time, flux, guessDict, freeParPlanet)
	# Readability of output data
	period = fitT.fitresultplanets['pnum0']['period']
	epoch = fitT.fitresultplanets['pnum0']['T0']
	k = fitT.fitresultplanets['pnum0']['rprs']
	rho = fitT.fitresultstellar['rho']
	duration = computeTransitDuration(period, rho, k)
	if not duration:
		duration = S.duration * 24
	# Calculating transit depth significance
	## fitT.transitmodel sometimes has a NaN value
	sigma = computePointSigma(time, flux, fitT.transitmodel, period, epoch, duration)
	depth = k ** 2
	significance = depth / sigma
	phase = getPhase(time, flux, period, epoch)[0]
	nTransitPoints = np.sum((-duration * 0.5 < phase) & (phase < duration * 0.5))
	SNR = significance * nTransitPoints**0.5
	return SNR, period, epoch, duration, depth, fitT.transitmodel, S.f_1, S.convolved_bls

def computePointSigma(time, flux, transitModel, period, epoch, duration):
	t2, f2 = removeTransits(time, flux, period, epoch, duration)
	mt2, mf2 = removeTransits(time, transitModel, period, epoch, duration)
	return np.nanstd(f2 - mf2)

def removeTransits(time, flux, period, epoch, duration):
	halfDur = 0.5 * duration / 24.
	bad = np.where(time < epoch - period + halfDur)[0]
	for p in np.arange(epoch, time[-1] + period, period):
		bad = np.append(bad, np.where((p - halfDur < time) & (time < p + halfDur))[0])
	good = np.setxor1d(range(len(time)), bad)
	return time[good], flux[good]

def computeTransitDuration(period, rho, k):
	b = 0.1							# Impact parameter (default value in ktransit)
	G = 6.67384e-11					# Gravitational constant
	P = period * 86400				# Period in seconds
	stellarDensity = rho * 1000
	rStarOverA = ((4 * np.pi**2) / (G * stellarDensity * P**2))**(1./3.)
	cosI = b * rStarOverA
	sinI = np.sqrt(1 - cosI**2)
	coeff = rStarOverA * np.sqrt((1+k)**2 - b**2) / sinI
	if coeff > 1:
		return 0
	else:
		duration = (P / np.pi) * np.arcsin(coeff)
		return duration / 3600		# Duration in hours


def findSecondary(time, flux, period, epoch, duration):
	t2, f2 = removeTransits(time, flux, period, epoch, duration)
	minp, maxp = period - 0.1, period + 0.1
	if t2[-1] - t2[0] == 0 or 1./maxp < 1./(t2[-1] - t2[0]):
		return (np.nan,)*5
	if minp < 0.5:
		minp = 0.5
	planetInfo = do_bls_and_fit(t2, f2, minp, maxp)
	return (t2,) + planetInfo[0:4] + (planetInfo[5],)

def computeOddEvenModels(time, flux, per, epo):
	gdOdd = {'period':	per * 2,
			 'T0':		epo}
	gdEven = {'period':	per * 2,
			  'T0':		epo + per}
	freeParPlanet = ['rprs']
	fitT_odd = fitModel(time, flux, gdOdd, freeParPlanet)
	fitT_even = fitModel(time, flux, gdEven, freeParPlanet)
	return fitT_odd, fitT_even

def main(filename):
	"""Fit a transit model to a lightcurve.
	1. Remove outliers.
	2. Detrend the data with a binwidth of 26 cadences. Since MAX_DURATION_HOURS = 6,
	   and 6 hours = ~13 cadences (ceiling of 12.245), detrending with a binwidth of double
	   this value will preserve all events with a duration of 13 cadences or less.
	3. Create an "S" object. (???) [1 is added to the flux to avoid a division by zero error]
	4. Run the BLS algorithm and Tom's transit fitting algorithm. Since the BLS can lock
	   on to an incorrect, shorter period event, I run it on four different minimum periods,
	   chosen somewhat arbitrarily. These results go into a dictionary sorted by the
	   calculated SNR of each fit, and the parameters which give the maximum SNR are used.
	5. Plot the original lightcurve, the BLS statistics from its minimum to maximum period,
	   and a phased lightcurve.
	6. Save the plot, and return a string containing the parameters of the fit.
	"""
	name = filename[-13:-4]
	time, flux = np.genfromtxt(filename, unpack = True)
	if np.all(np.isnan(flux)):
		return '%s\t\t%-8.6g\t%-8.6g\t%-8.6g\t%-4.3f\t%-8.6g\t%-8.6g' %((name,)+(np.nan,)*6)
	time, flux = outlierRemoval(time, flux)
	flux = medianDetrend(flux, 26)
	# Main transit search
	minPeriod = 0.5		# Limitations of BLS Fortran code
	maxPeriod = (time[-1] - time[0]) / 2.
	SNR, period, epoch, duration, depth, transitModel, period_guesses, \
	convolved_bls = do_bls_and_fit(time, flux, minPeriod, maxPeriod)
	# For the phase curves
	phase, phasedFlux = getPhase(time, flux, period, epoch)
	phaseModel, phasedFluxModel = getPhase(time, transitModel, period, epoch)
	# Secondary search
	secTime, secSNR, secPer, secEpoch, secDur, secModel = findSecondary(time, flux, period, epoch, duration)
	if secSNR > 5 and abs(period - secPer) < 0.05:
		secPhase, secPhaseModel = getPhase(secTime, secModel, secPer, epoch)
		idx = len(secPhase[secPhase < 0])
	else:
		secPhase, secPhaseModel, idx = [], [], 1
	# Odd/Even plot	
	fitT_odd, fitT_even = computeOddEvenModels(time, flux, period, epoch)
	phaseModel_odd, phasedFluxModel_odd = getPhase(time, fitT_odd.transitmodel, period * 2, epoch)
	phaseModel_even, phasedFluxModel_even = getPhase(time, fitT_even.transitmodel, period * 2, epoch + period)
	depthOdd = fitT_odd.fitresultplanets['pnum0']['rprs'] ** 2
	depthEven = fitT_even.fitresultplanets['pnum0']['rprs'] ** 2
	phaseOdd, fluxOdd = getPhase(time, flux, period * 2, epoch)
	phaseEven, fluxEven = getPhase(time, flux, period * 2, epoch + period)
	x1, x2 = -duration, duration
	y1, y2 = -3*np.std(fluxOdd), 3*np.std(fluxOdd)
	if min(fluxOdd) < y1:
		y1 = min(fluxOdd) - np.std(fluxOdd)
	# sigma = abs(depth1 - depth2) / sqrt(u1^2 + u2^2)
	durOdd = computeTransitDuration(period, fitT_odd.fitresultstellar['rho'], fitT_odd.fitresultplanets['pnum0']['rprs'])
	durEven = computeTransitDuration(period, fitT_odd.fitresultstellar['rho'], fitT_even.fitresultplanets['pnum0']['rprs'])
	sigma = computePointSigma(time, flux, transitModel, period, epoch, duration)
	nOddPoints = np.sum((-durOdd*0.5 < phaseOdd) & (phaseOdd < durOdd * 0.5))
	nEvenPoints = np.sum((-durEven*0.5 < phaseEven) & (phaseEven < durEven * 0.5))
	uOdd, uEven = sigma / np.sqrt(nOddPoints), sigma / np.sqrt(nEvenPoints)
	depthDiffSigma = abs(depthOdd - depthEven) / np.sqrt(uOdd**2 + uEven**2)
	if doPlot:
		gs = gridspec.GridSpec(3,2)
		ax1 = plt.subplot(gs[0,:])
		axOdd = plt.subplot(gs[1,0])
		axEven = plt.subplot(gs[1,1])
		ax3 = plt.subplot(gs[2,:])
		gs.update(wspace = 0, hspace = 0.5)
		ax1.plot(time, flux, 'k')
		y1, y2 = ax1.get_ylim()
		ax1.vlines(np.arange(epoch, time[-1], period), y1, y2, 
				   color = 'r', linestyles = 'dashed', linewidth = 0.5)
		ax1.axis([time[0], time[-1], y1, y2])
		ax1.set_title('kplr%s;    best period = %8.6g days;    SNR = %8.6g' %(name, period, SNR))
		ax1.set_xlabel('days')
		axOdd.set_ylabel('flux')
		axOdd.scatter(phaseOdd, fluxOdd, marker = '.', s = 1, color = 'k', alpha = 1)
		axOdd.plot(phaseModel_odd, phasedFluxModel_odd, 'r')
		axOdd.axhline(-depthOdd, x1, x2)
		axOdd.axis([x1,x2,y1,y2])
		axOdd.set_title('odd')
		axEven.scatter(phaseEven, fluxEven, marker = '.', s = 1, color = 'k', alpha = 1)
		axEven.plot(phaseModel_even, phasedFluxModel_even, 'r')
		axEven.axhline(-depthEven, x1, x2)
		axEven.yaxis.tick_right()
		axEven.axis([x1,x2,y1,y2])
		axEven.set_title('even')
		if secondary:
			plt.plot(secPhase[:idx], secPhaseModel[:idx], 'c')
			plt.plot(secPhase[idx:], secPhaseModel[idx:], 'c')
		ax3.scatter(phase, phasedFlux, marker = '.', s = 1, color = 'k')
		ax3.plot(phaseModel, phasedFluxModel, 'r')
		y1, y2 = -3*np.std(phasedFlux), 3*np.std(phasedFlux)
		if min(phasedFlux) < y1:
			y1 = min(phasedFlux) - np.std(phasedFlux)
		ax3.axis([phase[0], phase[-1], y1, y2])
		ax3.set_xlabel('phase [hours]')
		ax3.text(0.5, 1.25, 'depth diff sigma = %.3f' %depthDiffSigma, horizontalalignment = 'center',
			verticalalignment = 'center', transform = ax3.transAxes)
		if plotOption == 'save':
			plt.savefig(figureSaveLocation + '%s.png' %name, dpi = 200)
			plt.close()
		elif plotOption == 'show':
			plt.show()
	
	successString = '%s\t\t%-8.6g\t%-8.6g\t%-8.6g\t%-4.3f\t%-8.6g\t%-8.6g' \
					%(name, SNR, period, depth, epoch, duration, secSNR)
	return successString
	
def getResults():
	rfn = '/Users/Yash/Desktop/NASA/Summer2014/k2/changedWhichpix/run1/results.txt'
	names, periods = np.genfromtxt(rfn, usecols = (0,2), unpack = True)
	return names, periods
if __name__ == '__main__':
	# files = np.array(glob.glob('/Users/Yash/Desktop/NASA/Summer2014/k2/changedWhichpix/dataSVD/*.txt'))
	files = np.genfromtxt('/Users/Yash/Desktop/f0pcs.txt', dtype = 'str')
	title = '\t'.join(['name\t\t', 'SNR\t', 'period[days]', 'depth\t', 
					   'epoch[day]', 'duration[hours]', 'secondary SNR'])
	print(title)
	# Multiprocesses the code into a Pool of 7, while writing the results to 
	# resultsFilename as each iteration of the code completes. Also prints 
	# results to the console and gives an ETA.
	
	#-------------------------- MULTIPROCESSING --------------------------#
	# with open(resultsFilename, 'w') as rf:
	# 	rf.write(title + '\n')
	# p = Pool(7)
	# start = timeit.default_timer(); progressString = ''
	# for i, res in enumerate(p.imap_unordered(main, files), 1):
	# 	with open(resultsFilename, 'a') as rf:
	# 		rf.write(res + '\n')
	# 	avg = (timeit.default_timer() - start)/i
	# 	eta = (len(files) - i) * avg
	# 	sys.stdout.write('\b \b'*len(progressString))
	# 	print(res)
	# 	progressString = '%i/%i done, avg %3.2f sec per target, eta: %s' %(i, len(files), 
	# 		avg, pythonTime.strftime('%H:%M:%S', pythonTime.gmtime(eta)))
	# 	sys.stdout.write(progressString); sys.stdout.flush()
	# p.close()
	# p.join()
	# total = timeit.default_timer() - start
	# print('\ntotal elapsed time: %s' %pythonTime.strftime('%H:%M:%S', pythonTime.gmtime(total)))
