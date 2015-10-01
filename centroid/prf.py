
import numpy as np
import astropy.io.fits as pyfits
import os

__version__ = "$Id: prf.py 1964 2015-02-27 18:33:11Z fmullall $"
__URL__ = "$URL: svn+ssh://flux/home/fmullall/svn/kepler/py/prf.py $"


class KeplerPrf():
    """Return the expected PRF for a point source in the Kepler field
    based on mod out, and centroid.

    To get the recorded prf, use
    getPrfAtColRow()

    Other functions are in the works to map that PRF onto
    a given mask.

    Note:
    --------
    For speed, this function caches some data for each module to
    save reading it from file for each call. If you request the PRF
    from many mod-outs you may find your memory increases dramatically.

    """

    def __init__(self, path):
        self.path = path
        self.cache = dict()


    def getPrfForImage(self, img, hdr, mod, out, col, row):
        """Figure out bbox then call getPrfForBbox"""
        raise NotImplementedError("")
        #DO I need to worry about 1 based/ zero based here?
        c0 = float(hdr['1CRV4P'])
        r0 = float(hdr['2CRV4P'])

        img = self.getPrf(mod, out, img.shape, [c0, r0], \
            centroidOffsetColRow)
        return img


    def getPrfForBbox(self, mod, out, col, row, bboxIn):
        """Get the prf for an image described by a bounding box

        Input:
        -----------
        mod, out
            (ints) Which channel is a tgt on
        col, row
            (floats) Centroid for which to generate prf
        bboxIn
            (int array). Size of image to return. bbox
            is a 4 elt array of [col0, col1, row0, row1]

        Returns:
        ----------
        A 2d numpy array of the computed prf at the
        given position.

        Notes:
        ------------
        If bbox is larger than the prf postage stamp returned,
        the missing values will be filled with zeros. If
        the bbox is smaller than the postage stamp, only the
        requestd pixels will be returned
        """

        bbox = np.array(bboxIn).astype(int)
        nColOut = bbox[1] - bbox[0]
        nRowOut = bbox[3] - bbox[2]
        imgOut = np.zeros( (nRowOut, nColOut) )

        #Location of origin of bbox relative to col,row.
        #This is usually zero, but need not be.
        colOffsetOut = (bbox[0] - np.floor(col)).astype(np.int)
        rowOffsetOut = (bbox[2] - np.floor(row)).astype(np.int)
        #assert(colOffsetOut < 0)
        #assert(rowOffsetOut < 0)

        interpPrf = self.getPrfAtColRow(mod, out, col, row)
        nRowPrf, nColPrf = interpPrf.shape
        colOffsetPrf = -np.floor(nColPrf/2.).astype(np.int)
        rowOffsetPrf = -np.floor(nRowPrf/2.).astype(np.int)

        di = colOffsetPrf - colOffsetOut
        i0 = max(0, -di)
        i1 = min(nColOut-di , nColPrf)
        if i1 <= i0:
            raise ValueError("Central pixel column not in bounding box")
        i = np.arange(i0, i1)
        assert(np.min(i) >= 0)

        dj = rowOffsetPrf - rowOffsetOut
        j0 = max(0, -dj)
        j1 = min(nRowOut-dj, nRowPrf)
        if j1 <= j0:
            raise ValueError("Central pixel row not in bounding box")
        j = np.arange(j0, j1)
        assert(np.min(j) >= 0)

        #@TODO: figure out how to do this in one step
        for r in j:
            imgOut[r+dj, i+di] = interpPrf[r, i]

        #imgOut[i+di, j+dj] = interpPrf[i,j]
        return imgOut


    def old_getPrfAtColRow(self, mod, out, col, row):
        """
        Computes the prf without using the cache. Left behind for
        debugging purposes.
        """
        fullPrfArray = self.getSubSampledPrfs(mod, out)
        regPrfArray = self.getRegularlySampledPrfs(fullPrfArray, col,row)
        bestRegPrf = self.interpolateRegularlySampledPrf(regPrfArray, col, row)

        return bestRegPrf


    def getPrfAtColRow(self, mod, out, col, row):
        """Compute the model prf for a given module, output, column row

        This is the workhorse function of the class. For a given mod/out,
        loads the subsampled PRFs at each corner of the mod-out, extracts
        the appropriate image for the given subpixel position, then interpolates
        those 4 images to the requested column and row.
        """

        #Load subsampled PRFs from cache, or from file if not previously read in
        key = "%02i-%02i" %(mod, out)
        if key not in self.cache:
            self.cache[key] = self.getSubSampledPrfs(mod, out)

        fullPrfArray = self.cache[key]
        regPrfArray = self.getRegularlySampledPrfs(fullPrfArray, col,row)
        bestRegPrf = self.interpolateRegularlySampledPrf(regPrfArray, \
            col, row)

        return bestRegPrf


    def getSubSampledPrfs(self, mod, out):
        fullPrfArray = []
        for i in range(1,6):
            tmp = self.readPrfFile(mod, out, i)
            fullPrfArray.append(tmp)
        return fullPrfArray


    def readPrfFile(self, mod, out, position):
        """

        position    (int) [1..5]
        """

        filename = "kplr%02i.%1i_2011265_prf.fits" %(mod, out)
        filename = os.path.join(self.path, filename)

        img = pyfits.getdata(filename, ext=position)
        return img


    def getRegularlySampledPrfs(self, fullPrfArray, col, row):
        regArr = []
        for i in range(5):
            tmp = self.getSingleRegularlySampledPrf(fullPrfArray[i], \
                col, row)
            regArr.append(tmp)
        return regArr


    def getSingleRegularlySampledPrf(self, singleFullPrf, col, row):
        """Extract out an 11x11* PRF for a given column and row for
        a single 550x550 representation of the PRF

        Doesn't interpolate across prfs just takes
        care of the intrapixel variation stuff

        Steve stored the prf in a funny way. 50 different samplings
        of an 11x11 pixel prf*, each separated by .02 pixels are
        stored in a 550x550 grid. Elements in the grid 50 spaces
        apart each refer to the same prf.


        *Sometimes a 15x15 pixel grid

        Notes:
        ---------
        None of the details of how to write this function are availabe
        in the external documents. It was all figured out by trial and error.
        """
        gridSize = 50.   #Defined in documentation

        #The sub-pixel image from (0.00, 0.00) is stored at x[49,49].
        #Go figure.
        colIndex = ((1-np.remainder(col, 1)) * gridSize).astype(np.int32)-1
        rowIndex = ((1-np.remainder(row, 1)) * gridSize).astype(np.int32)-1

        nColOut, nRowOut = singleFullPrf.shape
        nColOut /= float(gridSize)
        nRowOut /= float(gridSize)

        iCol = colIndex + (np.arange(nColOut)*gridSize).astype(np.int)
        iRow = rowIndex + (np.arange(nRowOut)*gridSize).astype(np.int)

        #Don't understand why this must be a twoliner
        tmp = singleFullPrf[iRow, :]
        return tmp[:,iCol]


    def interpolateRegularlySampledPrf(self, regPrfArray, col, row):

        #See page 2 of PRF_Description_MAST.pdf for the storage
        #order
        p11, p12, p21, p22, pMid = regPrfArray

        #See instrument hand bok, page 49, Section 4.5
        nCol, nRow = 1099, 1023
        col = int(col) / float(nCol)
        row = int(row) / float(nRow)

        #Intpolate across the rows
        tmp1 = p11 + (p21 - p11) * col
        tmp2 = p12 + (p22 - p12) * col

        #Interpolate across the columns
        out = tmp1 + (tmp2-tmp1)*row
        return out




    def mapPrfToImg(self, bestRegPrf, imgSizeRC, imgOffsetCR, \
            centroidCR):
        """Place a rectangular apeture over the prf

        Note this will fail if the image aperture doesn't wholly
        overlap with bestRegPrf
        """

        #Location of brightest pixel in PRF img. Usually 5,5,
        #but sometimes 7,7
        midPrfRow, midPrfCol = np.floor(np.array(bestRegPrf.shape) / 2.)
        nRowImg, nColImg = imgSizeRC

        #Map midPrf so it lies on the same pixel within the image
        #centroidCR does
        xy = np.array(centroidCR) - np.array(imgOffsetCR)
        #c1 = midPrfCol-np.floor(centroidCR[0]) + 1
        #c1 = np.arange(c1, c1+imgSizeRC[1], dtype=np.int32)

        #r1 = midPrfRow-np.floor(centroidCR[1]) + 1
        #r1 = np.arange(r1, r1+imgSizeRC[0], dtype=np.int32)
        c1 = midPrfCol-xy[0] + 1
        c1 = np.arange(c1, c1+imgSizeRC[1], dtype=np.int32)

        r1 = midPrfRow-xy[1] + 1
        r1 = np.arange(r1, r1+imgSizeRC[0], dtype=np.int32)

        tmp = bestRegPrf[r1, :]
        return tmp[:,c1]

