import numpy as np
import pylab as plt
import time

def int_nonlin(lim, c):
    norms = [5.0, 3.0, 7.0, 2.0]
    facs = [4.0, 2.0, 4.0]
    sqlim = lim * lim
    tmplim = np.sqrt(1.0 - sqlim)
    tmplimquart = np.sqrt(tmplim)
    tmplim3quart = tmplimquart * tmplimquart * tmplimquart
    workf = lambda x, y, z: (x + (y * (z - x))) / z
    return sqlim - (c[0]*workf(facs[0]*tmplimquart, sqlim, norms[0] ) + \
                    c[1]*workf(facs[1]*tmplim, sqlim, norms[1]) + \
                    c[2]*workf(facs[2]*tmplim3quart, sqlim, norms[2]) + \
                    c[3]*sqlim*sqlim/norms[3])
    
def int_quad(lim, c):
    sqlim = lim * lim
    tmplim = 2.0 * np.sqrt(1.0 - sqlim)
    return sqlim - (c[0]*(tmplim + (sqlim * (3.0 - tmplim)))/3.0  + \
                    c[1]*sqlim*sqlim*0.5)

    
def transit_model_small_planet_nonlinear(z, k, c1, c2, c3, c4):
    """ Model a transit light curve in the small planet approximation
        of Mandel & Agol using the nonlinear limb darkening law
        Works up to k=Rp/Rstar =0.5, but recommend much smaller values
        AUTHOR: Christopher J. Burke
            based upon the mathematical description of Mandel & Agol
    """
    lc = np.full_like(z, 1.0)
    c0 = 1.0 - c1 - c2 - c3 - c4
    omega = c0/4.0 + c1/5.0 + c2/6.0 + c3/7.0 + c4/8.0
    # check if one should use quad instead of nonlinear
    if c1 == 0.0 and c3 == 0.0:
        cs = [c2, c4]
        intuse = int_quad
        intzero = -c2 * 2.0 / 3.0
        intone = 1.0 - c2 - 0.5*c4
    else:
        cs = [c1, c2, c3, c4]
        intuse = int_nonlin
        intzero = -(c1*4.0/5.0 + c2*2.0/3.0 + c3*4.0/7.0)
        intone = 1.0 - c1 - c2 - c3 - 0.5*c4
    # Do cases where k<=z<=1-k
    idx = np.where(np.logical_and(z <= 1.0-k, z >= k))[0]
    if idx.size > 0:
        zz = z[idx]
        if k.size > 1:
            kk = k[idx]
        else:
            kk = k
        iStar = (intuse(zz+kk, cs) - \
                 intuse(zz-kk, cs)) / \
                 (4.0 * kk * zz)
        lclc = 1.0 - kk * kk * iStar / 4.0 / omega
        lc[idx] = lclc
    # Do cases where 1-k<z<=1+k
    idx = np.where(np.logical_and(z < 1.0+k, z > 1.0-k))[0]
    if idx.size > 0:
        zz = z[idx]
        if k.size > 1:
            kk = k[idx]
        else:
            kk = k
        iStar = (intone - \
                 intuse(zz-kk, cs)) / \
                (1.0 - (zz-kk)*(zz-kk))
        lclc = 1.0 - iStar / 4.0 / omega / np.pi * (kk**2 * \
               np.arccos((zz - 1.0)/kk) - (zz - 1.0)* \
               np.sqrt(kk**2 - (zz-1.0)**2))
        lc[idx] = lclc
    # Do case where 0<=z<k
    idx = np.where(z < k)[0]
    if idx.size > 0:
        zz = z[idx]
        if k.size > 1:
            kk = k[idx]
        else:
            kk = k
        lk = zz + kk
        iStarl = (intuse(lk,  cs) - \
                  intzero) 
        rk = kk - zz
        iStarr = (intuse(rk,  cs) - \
                  intzero) 
        lclc = 1.0 - kk * kk * (iStarl + iStarr) / (lk*lk + rk*rk) / 4.0 / omega
        lc[idx] = lclc
    return lc

def transit_model_small_planet_quad(z, k, u1, u2):
    """ Model a transit light curve in the small planet approximation
        of Mandel & Agol using the quadratic limb darkening law
        Works up to k=Rp/Rstar =0.5, but recommend much smaller values
        AUTHOR: Christopher J. Burke
            based upon the mathematical description of Mandel & Agol

    """
    lc = transit_model_small_planet_nonlinear(z, k, 0.0, u1+2.0*u2, 0.0, -u2) 
    return lc


# Run the test of a small planet in gaussian noise
if __name__ == "__main__":

    usek = np.array([0.1])
    z = np.linspace(0.0,1.0+usek,1000000)
    t0 = time.time()
    lc = transit_model_small_planet_quad(z, usek, 0.40, 0.27)
    print time.time() - t0, "seconds wall time"
    t0 = time.time()
    lc2 = transit_model_small_planet_nonlinear(z, usek, 0.432, 0.1787, 0.3969, -0.2577)
    print time.time() - t0, "seconds wall time"
    plt.plot(z,lc,'.b')
    plt.plot(z,lc2,'.r')
    plt.show()


