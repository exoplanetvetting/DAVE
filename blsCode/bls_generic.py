from __future__ import division, print_function
import numpy as np
import yash_bls 

def doSearch(time, flux, minPeriod, maxPeriod, ):

    SNR, period, epoch, duration, depth, transitModel, period_guesses, \
            convolved_bls = yash_bls.do_bls_and_fit(time, flux, minPeriod, maxPeriod)

    phase, phasedFlux = yash_bls.getPhase(time, flux, period, epoch)
    phaseModel, phasedFluxModel = yash_bls.getPhase(time, transitModel, period, epoch)

    secTime, secSNR, secPer, secEpoch, secDur, secModel = yash_bls.findSecondary(time, flux, period, epoch, duration)
    if secSNR > 5 and abs(period - secPer) < 0.05:
        secPhase, secPhaseModel = yash_bls.getPhase(secTime, secModel, secPer, epoch)
        idx = len(secPhase[secPhase < 0])
    else:
        secPhase, secPhaseModel, idx = [], [], 1


    # Odd/Even plot 
    fitT_odd, fitT_even = yash_bls.computeOddEvenModels(time, flux, period, epoch)
    phaseModel_odd, phasedFluxModel_odd = yash_bls.getPhase(time, fitT_odd.transitmodel, period * 2, epoch)
    phaseModel_even, phasedFluxModel_even = yash_bls.getPhase(time, fitT_even.transitmodel, period * 2, epoch + period)
    depthOdd = fitT_odd.fitresultplanets['pnum0']['rprs'] ** 2
    depthEven = fitT_even.fitresultplanets['pnum0']['rprs'] ** 2
    phaseOdd, fluxOdd = yash_bls.getPhase(time, flux, period * 2, epoch)
    phaseEven, fluxEven = yash_bls.getPhase(time, flux, period * 2, epoch + period)
    x1, x2 = -duration, duration
    y1, y2 = -3*np.std(fluxOdd), 3*np.std(fluxOdd)
    if min(fluxOdd) < y1:
        y1 = min(fluxOdd) - np.std(fluxOdd)
    # sigma = abs(depth1 - depth2) / sqrt(u1^2 + u2^2)
    durOdd = yash_bls.computeTransitDuration(period, fitT_odd.fitresultstellar['rho'], fitT_odd.fitresultplanets['pnum0']['rprs'])
    durEven = yash_bls.computeTransitDuration(period, fitT_odd.fitresultstellar['rho'], fitT_even.fitresultplanets['pnum0']['rprs'])
    sigma = yash_bls.computePointSigma(time, flux, transitModel, period, epoch, duration)
    nOddPoints = np.sum((-durOdd*0.5 < phaseOdd) & (phaseOdd < durOdd * 0.5))
    nEvenPoints = np.sum((-durEven*0.5 < phaseEven) & (phaseEven < durEven * 0.5))
    uOdd, uEven = sigma / np.sqrt(nOddPoints), sigma / np.sqrt(nEvenPoints)
    depthDiffSigma = abs(depthOdd - depthEven) / np.sqrt(uOdd**2 + uEven**2)

    return locals


def plotSearch():
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




if __name__ == '__main__':

    filename = 'wasp47.csv'

    time, flux = np.genfromtxt(filename, unpack = True, delimiter=',')

    mask = np.isfinite(flux) * np.isfinite(time)
    time = time[mask]
    flux = flux[mask]


    time, flux = yash_bls.outlierRemoval(time, flux)
    flux = yash_bls.medianDetrend(flux, 26)

# Main transit search
    minPeriod = 0.5     # Limitations of BLS Fortran code
    maxPeriod = (time[-1] - time[0]) / 3.

    outs = doSearch(time, flux, minPeriod, maxPeriod, )


    
    #successString = '%s\t\t%-8.6g\t%-8.6g\t%-8.6g\t%-4.3f\t%-8.6g\t%-8.6g' \
     #               %(name, SNR, period, depth, epoch, duration, secSNR)