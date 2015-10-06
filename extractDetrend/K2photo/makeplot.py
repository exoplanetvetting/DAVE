from __future__ import division, print_function

import numpy as np
import matplotlib.pyplot as plt

import pyfits
import glob

import extract_lc
import progressbar

params = {'backend': 'png',
            'axes.linewidth': 2.5,
            'axes.labelsize': 24,
            'axes.font': 'sans-serif',
            'axes.fontweight' : 'bold',
            'text.fontsize': 12,
            'legend.fontsize': 16,
            'xtick.labelsize': 16,
            'ytick.labelsize': 16,
#            'text.usetex': True,
#            'font.family': 'sans-serif',
            'font.sans-serif': 'Helvetica',
            'ps.useafm': True,
            'pdf.use14corefonts': True,
            'ps.fonttype': 42,
            'legend.markersize': 200}
plt.rcParams.update(params)

def rolling_window(a, window):
    shape = a.shape[:-1] + (a.shape[-1] - window + 1, window)
    strides = a.strides + (a.strides[-1],)
    return np.lib.stride_tricks.as_strided(a, shape=shape, strides=strides)


if __name__ == '__main__':
    dirname = '/Users/tom/Projects/Debra_data/data/lcfiles/'
    files = glob.glob(dirname + '*.txt')

    cdpparr = np.zeros(len(files))

    bar = progressbar.ProgressBar(maxval=len(files), \
    widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage(),
            progressbar.ETA()])

    bar.start()


    for i,f in enumerate(files):

        bar.update(i+1)

        fig, (ax1,ax2,ax3) = plt.subplots(3,1, sharex=True, sharey=True,
            figsize=[9,12])

        t = np.genfromtxt(f).T

        ax1.plot(t[0],t[1],c='r')
        ax2.plot(t[0],t[2],c='b')
        ax3.plot(t[0],t[3],c='g')

        #calc CDPP
        cdpp = 1.E6 * np.median(
            np.std(rolling_window(t[2], 13) / np.sqrt(13), 1))
        cdpparr[i] = cdpp

        ax1.set_xlim([45.5,77.0])

        ax1.set_title('Decor lc,')
        ax2.set_title('Decor + medfilt lc, CDPP = {}'.format(cdpp))
        ax3.set_title('Decor signal')


        plt.tight_layout()

        savename = f.split('/')[-1].split('.')[0]
        plt.savefig('{}../figs/{}.png'.format(dirname,savename))

        plt.close('all')

    bar.finish()
    np.savetxt('cdpp.txt',np.array([files,cdpparr],dtype=None).T,fmt='%s')














