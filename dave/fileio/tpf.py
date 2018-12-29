
#Tools for manipulating target pixel files
#$Id: tpf.py 1305 2013-03-03 21:51:25Z fmullall $
#$URL: svn+ssh://flux/home/fmullall/svn/kepler/py/tpf.py $
#Fergal Mullally


import numpy as np


def getTargetPixelArrayFromFits(fits, hdr, column='FLUX'):
    nRows = int(hdr['NAXIS2'])
    #All imgs are the same size, so it doesn't matter which TDIM
    #I pick here.
    shape = eval(hdr['TDIM5'])
    #The reverse puts the columns on the x-axis and row
    #on the y-axis when you call mp.imshow()
    shape = tuple(reversed(shape))

    tpfArray = np.empty((nRows, shape[0], shape[1]))

    for i, cadence in enumerate(fits[column]):
        tpfArray[i,:,:] = cadence.reshape(shape)

    return tpfArray


def lightCurveFromTpfCube(cube, mask=None):
    """Extract a lightcurve from a TPF cube.

    Inputs:
    cube    (3d numpy array) As returned by getTargetPixelArrayFromFits()
    mask    (2d numpy array) Which pixels to use in extraction. Default
            is all finite pixels

    Returns
    A 1d numpy array

    Notes:
    Pixels marked Nan (i.e outside the aperture) are ignored regardless
    of the mask value.

    Default mask includes all finite pixels.

    To extract only the optimal aperture use
    mask = np.where(mask == 3, 1, 0)
    lightCurveFromTpfCube(cube, mask)
    """

    if mask is None:
        mask = np.ones( cube[0,:,:].shape)

    n, nR, nC = cube.shape
    out = np.zeros(n)
    for i in range(n):
        img = cube[i, :, :]
        idx = np.bitwise_and(mask > 0, np.isfinite(img))
        out[i] = np.sum(img[idx])

    return out

