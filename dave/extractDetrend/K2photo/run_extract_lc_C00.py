from __future__ import division, print_function

import numpy as np
import matplotlib.pyplot as plt

import pyfits
import glob

import extract_lc
import progressbar

if __name__ == '__main__':
    filenames = glob.glob(
        '/Volumes/K2_ENG/C0_all/ktwo202*lpd-targ.fits.gz')

    bar = progressbar.ProgressBar(maxval=len(filenames), \
    widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage(),
            progressbar.ETA()])

    bar.start()

    for i,filename in enumerate(filenames):
        #print('{} of {}'.format(i,len(filenames)))
        bar.update(i+1)

        time,lc,xbar,ybar,regnum = extract_lc.run_C0_data_extract(
            filename, bg_cut=4)

        if regnum == 0:
            continue

        #mask = np.ones(len(lc),dtype=bool)
        #mask[range(11,len(lc),12)] = False
        #time,lc,xbar,ybar = time[mask], lc[mask], xbar[mask], ybar[mask]


        #lets flatten the light curve
        flatlc = extract_lc.medfilt(time,lc,window=1.)
        zpt = len(time)%300.

        try:
            outflux, correction, thr_cad = extract_lc.run_C0_detrend(
                time,flatlc,xbar,ybar,cadstep=300)
        except ValueError:
            print('{} had an ValueError'.format(filename))
            continue

        not_thr = ~thr_cad

        corflux = (lc[zpt:][not_thr]/
                    np.median(lc[zpt:][not_thr])/
                    correction[not_thr])
        corflatflux = (flatlc[zpt:][not_thr]/
                    np.median(flatlc[zpt:][not_thr])/
                    correction[not_thr])


        outname = filename.split(
            '/')[-1].split(
            'pd-targ.fits.gz')[0]
        np.savetxt(
            '/Users/tom/Projects/Debra_data/data/{}.txt'.format(
                outname),
            np.array([time[zpt:][not_thr],
                        corflux,corflatflux,
                        correction[not_thr]]).T)


        #plt.figure()
        #plt.plot(time[zpt:][not_thr],
        #    lc[zpt:][not_thr]/np.median(lc[zpt:][not_thr])
        #    /correction[not_thr])
        #plt.scatter(time[zpt:],
        #    flatlc[zpt:],s=2,color='b')


    bar.finish()




