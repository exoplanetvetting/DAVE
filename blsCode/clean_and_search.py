# import sys
# sys.path.append('/Users/Yash/svn_code/tom_code')
# import medium_filter
# import sgfilt
import numpy as np
# import untrendy
# from pybls import bls
import scipy.ndimage.filters
import matplotlib.pyplot as plt
import bls
# import filter
# import kplr
from numpy.random import choice
import logging
logging.disable('WARNING')

class Clean(object):
    def __init__(self):
        self.client = kplr.API()
    def get_data(self,kic,pdc=False):
        star = self.client.star(kic)
        lcs = star.get_light_curves(short_cadence=False)
        time, flux, ferr, quality,quarter = [], [], [], [], []
        for lc in lcs:
            try:
                with lc.open() as f:
                    # The lightcurve data are in the first FITS HDU.
                    hdu_data = f[1].data
                    time = np.r_[time,hdu_data["time"]]
                    # fluxval = hdu_data["sap_flux"]
                    # fluxmed = np.median(fluxval)
                    # flux = np.r_[flux,fluxval / fluxmed]
                    # ferr = np.r_[ferr,hdu_data["sap_flux_err"] / fluxmed]
                    if pdc:
                        flux = np.r_[flux,hdu_data["pdcsap_flux"]]
                        ferr = np.r_[ferr,hdu_data["pdcsap_flux_err"]]
                    else:
                        flux = np.r_[flux,hdu_data["sap_flux"]]
                        ferr = np.r_[ferr,hdu_data["sap_flux_err"]]
                    quality = np.r_[quality,hdu_data["sap_quality"]]
                    quarter = np.r_[quarter,f[0].header["QUARTER"] +
                        np.zeros(len(hdu_data["time"]))]
            except:
                # perhaps data is proprietary
                pass
        self.time = time
        self.flux = flux
        self.ferr = ferr
        self.quality = quality
        self.quarter = quarter
        self.byebyebaddata2()

    def gal_med(self,window=48):
        self.cflux = medium_filter.detrend_using_mf(self.flux,size=window)

    def byebyebaddata2(self):
        finite = np.isfinite(self.flux)
        qflag = self.quality == 0
        mask = np.logical_and(finite,qflag)
        self.flux[~mask] = np.nan
        self.ferr[~mask] = np.nan

    def read_data(self,time,flux,ferr,quality,quarter):
        self.time = time
        self.flux = flux
        self.ferr = ferr
        self.quality = quality
        self.quarter = quarter

    def remove_post_earth_point(self,time,flux,window=2.0):
        ep_idx = [i for i,x in enumerate(
            self.quality) if x in [1,2,8]]
        for tep in time[ep_idx]:
            trange = np.logical_and(time >= tep,time < tep + 2.0)
            flux[trange] = np.nan
        return flux

    def byebyebaddata(self):
        finite = np.isfinite(self.flux)
        qflag = self.quality == 0
        mask = np.logical_and(finite,qflag)
        self.time = self.time[mask]
        self.flux = self.flux[mask]
        self.ferr = self.ferr[mask]
        self.quality = self.quality[mask]
        self.quarter = self.quarter[mask]

    def norm_by_quarter(self):
        for i in np.unique(self.quarter):
            medold = np.median(
                self.flux[self.quarter == i])
            self.flux[self.quarter == i] /= medold
            self.ferr[self.quarter == i] /= medold

    def medfilt(self,window=2.0):
        """
        window in days
        """
        self.cflux = self.flux / untrendy.median(
            self.time,self.flux,dt=window)
        self.cferr = self.ferr

    def untrend(self):
        self.cflux, self.cferr = untrendy.untrend(
            self.time,self.flux,self.ferr)

    def MAD(self,xx):
        """Median Absolute Deviation
        """
        med=np.median(xx,0)
        absdev=np.abs(np.subtract(xx,med))
        mad=np.median(absdev,0)
        return 1.48 * mad

    def sig_clip_special(self,sigma=5.,maxiter=5):
        flux2 = np.copy(self.cflux)
        time2 = np.copy(self.time)
        ferr2 = np.copy(self.ferr)
        i = 0
        while i < maxiter:
            med_val = np.median(flux2)
            std_val = self.MAD(flux2) * sigma
            inc_data = np.logical_and(flux2 > (med_val-std_val),
                flux2 < (med_val+std_val))
            if len(flux2) == len(flux2[inc_data]):
                break
            flux2 = flux2[inc_data]
            time2 = time2[inc_data]
            err2 = err2[inc_data]
        self.ctime = time2
        self.cflux = flux2
        self.cferr = ferr2

class Search(object):
    def __init__(self,time,flux,ferr):
        self.time = time
        self.flux = flux
        self.ferr = ferr

    def do_bls(self,minp=200,maxp=500):
        """
        inputs to bls are:
        p,bper,bpow,depth,qtran,in1,in2 = eebls(
            t,x,e,nf,fmin,df,nb,qmi,qma,[n])
        where:
        t is time
        x is flux
        e is err2
        nf is number of frequencies (5000)
        fmin is minimum freq
        df is the stepsize ()
        nb is numbe rof bins
        qmin frac in transit
        qmax frac in transit
        compute_bls(time,
                    lc, df = freq_step,
                    nf =  (1./min_period)/freq_step, nb = 1400,
                    qmi = float(min_duration_hours)/24./450.,
                    qma = float(max_duration_hours)/24./300.,
                    fmin = (1./(float(max_period)*1.1)))
        """
        maxf = 1. / minp
        minf = 1. / maxp
        nf = 1
        df = 0
        samp = (15./1440.)
        perarr = np.arange(minp,maxp,samp)
        freqarr = 1. / perarr
        arrsize = len(freqarr)
        arr = np.zeros([arrsize,7])
        for i, (per,freq) in enumerate(zip(perarr,freqarr)):
            tt = (np.max(time)-np.min(time))
            tdur = (10 * (per / 365.256)**(1./3.)) / 24.
            qmin = qmax = tdur / tt
            nb = (tt/tdur) * 4
            p,bper,bpow,depth,qtran,in1,in2 = bls.bls.eebls(
            time,flux,ferr,nf,freq,df,nb,qmin,qmax)
            arr[i] = p,bper,bpow,depth,qtran,in1,in2
            if i%1000 == 0:
                print i
        self.arr= arr

    def compute_bls(self,time, lc, df, nf, nb, qmi, qma, fmin,norm=True):
        bls, f_1, nb = self.BLS(time, lc, df, nf, nb, qmi, qma, fmin)
        usethisbls = np.array(bls[0])
        if norm:
            boot_loop = 9
            boot_bls = np.empty([boot_loop,nf])
            lc_rand = choice(lc,size=len(time))
            for i in xrange(9):
                bls_r, f_1_r, nb_r = self.BLS(time,
                    lc_rand, df, nf, nb, qmi, qma, fmin)
                boot_bls[i] = bls_r[0]
            usethisbls = bls[0] / np.median(boot_bls,axis=0)
        bper = bls[1]
        bpow = bls[2]
        depth = bls[3]
        qtran = bls[4]
        duration = bper*qtran
        in1 = bls[5]
        in2 = bls[6]
        phase1 = in1/float(nb)
        phase2 = in2/float(nb)
        convolved_bls = scipy.ndimage.filters.gaussian_filter(
            usethisbls, 2.0)
        peak = np.r_[True, convolved_bls[1:] > convolved_bls[:-1]] & \
            np.r_[convolved_bls[:-1] > convolved_bls[1:], True]
        '''Sort peaks'''
        sel_peaks = np.sort(convolved_bls[peak])
        sel_peaks = sel_peaks[-1:]
        '''locate highest peak'''
        periods = f_1[np.where(convolved_bls == sel_peaks)]
        '''calculate number of transits, epoch, ingress and egress times'''
        t_number = int((max(time) - min(time)) / bper)
        epoch = time[0] + phase1*bper
        ingresses = np.zeros(t_number)
        egresses = np.zeros(t_number)
        for n in range(0,t_number):
            ingresses[n] = (epoch + bper*n) - 0.2
            egresses[n] = epoch + bper*n + duration + 0.2 # add a margin each side
        approx_duration = egresses[0] - ingresses[0]
        return ingresses, egresses, t_number, epoch, periods, bper, bpow, depth, qtran, \
            duration, f_1, convolved_bls, approx_duration

    def BLS(self,time, lc, df = 0.0001, nf = 500,  nb = 200, qmi = 0.01,\
        qma = 0.8, fmin = (1./(400.0*1.1))):
        diffs = time[1:] - time[:-1]
        u = np.ones(len(time))
        v = np.ones(len(time))
        BLS = bls.eebls(time, lc, u, v, nf, fmin, df, nb, qmi, qma)
        f = fmin + (np.arange(len(BLS[0])))*df
        return BLS, 1/f, nb

    def get_qf(self,time,flux,period,epoch):
        date1 = (time - epoch) + 0.5*period
        phi1 = (((date1 / period) - np.floor(date1/period)) * 24. * period) - 12*period
        q1 = np.sort(phi1)
        f1 = (flux[np.argsort(phi1)])
        return q1, f1

    def do_bls2(self,min_period=120.,max_period=500.,
        min_duration_hours=5,max_duration_hours=20,
        freq_step = 0.0000001,nbins = 1000,norm=True,
        doplot=True):
        finite = np.isfinite(self.flux)
        time = self.time[finite]
        lc = self.flux[finite] / np.median(self.flux[finite])
        # min_period = 120.
        # max_period = 500.
        # freq_step = 0.0000001
        # min_duration_hours = 5
        # max_duration_hours = 20
        # nbins = 1000
        ofac = 64.
        nyq = 24.
        freq1 = (1./max_period)
        freq2 = (1./min_period)
        freqdiff = freq2 - freq1
        steps=int(ofac*freqdiff*len(time) / nyq)
        df=freqdiff/steps
        ingresses, egresses, t_number, epoch, periods, bper, bpow, depth, qtran, \
        duration, f_1, convolved_bls, approx_duration = self.compute_bls(time, \
                lc, df = df, \
                nf =  steps, nb = nbins, \
                qmi = float(min_duration_hours)/24./max_period, \
                qma = float(max_duration_hours)/24./min_period, \
                fmin = freq1,norm=norm)
        if doplot:
            fig = plt.figure()
            plt.subplot(2,1,1)
            plt.plot(time,lc,".k")
            for i in range(len(ingresses)):
                plt.axvline(ingresses[i], color = 'c')
                plt.axvline(egresses[i], color = 'c')
            plt.ylabel('Flux')
            plt.xlabel('Time')
            plt.subplot(2,1,2)
            plt.plot(f_1, convolved_bls)
            plt.axvline(bper, color = 'r')
            plt.xlim(min(f_1), max(f_1))
            plt.xlabel('Period')
            for i in range(2, 10):
                plt.axvline(i*bper, color = 'y')
            plt.axvline(bper/2., color = 'y')
            plt.axvline(3*bper/2., color = 'y')
            fig2 = plt.figure()
            q1,f1 = self.get_qf(time,lc,bper,epoch)
            limit = np.logical_and(q1 > -25,q1 < 25)
            plt.scatter(q1[limit],f1[limit],s=1)
            plt.xlim([-25,25])
            for i in xrange(-25,25,2):
                bin = np.logical_and(q1 > i,q1 <= i+2)
                fval = f1[bin]
                fer = np.std(fval) / np.sqrt(len(fval))
                plt.errorbar(i+1,np.mean(fval),yerr=fer,
                    fmt='o',color='b',ms=5)
        # print 'period = ', bper
        self.ingresses = ingresses
        self.egresses = egresses
        self.t_number = t_number
        self.epoch = epoch
        self.periods = periods
        self.bper = bper
        self.bpow = bpow
        self.depth = depth
        self.qtran = qtran
        self.duration = duration
        self.f_1 = f_1
        self.convolved_bls = convolved_bls
        self.approx_duration = approx_duration
        
if __name__ == '__main__':
    #datafile = '/Users/tom/scratch/faketransits/ascii_inj/k2571748.dat'
    #time,flux,ferr,quality,quarter = np.loadtxt(
    #    datafile,unpack=True)
    obj = Clean()
    #obj.read_data(time,flux,ferr,quality,quarter)
    #obj.byebyebaddata()
    #obj.norm_by_quarter()
    #obj.medfilt()
    kic = 8478994
    obj.get_data(kic)
    #obj.gal_med(window=101)
    obj.norm_by_quarter()
    obj.medfilt()
    obj.cflux2 = obj.remove_post_earth_point(obj.time,obj.cflux)
    sea = Search(obj.time,(obj.cflux2),obj.ferr)
    sea.do_bls2()
    #p,bper,bpow,depth,qtran,in1,in2 = eebls(
    #        t,x,e,nf,fmin,df,nb,qmi,qma,[n])
    #t = bls.bls.eebls(obj.time,obj.cflux,obj.cferr,
    #    1./200./0.0000005,1./500.,0.0000005,4000.,5./24./500.,15./24./200.)
