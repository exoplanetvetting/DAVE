"""
Created on Sun Dec  2 14:12:33 2018

@author: fergal
"""

from __future__ import print_function
from __future__ import division

from dave.diffimg.AbstractPrfLookup import AbstractPrfLookup
from pdb import set_trace as debug
import scipy.io as spio
from glob import glob
import numpy as np
import os

npmap = lambda f, x: np.array(list(map(f, x)))



class TessPrf(AbstractPrfLookup):
    """Interpolate a TESS PRF image

    The TESS mission makes a model PRF available to the community through the MAST Archive.
    The model prf is evaluated at a grid of locations across each CCD, and for a number
    of subpixel positions at each grid location. This class contains the logic for
    extracting the PRF for an arbitrary location within a CCD.

    The two public methods of this class are
    .. code-block:: python
        TessPrf.getPrfAtColRow(col, row, ccd, camera, sector))
        TessPrf.getPrfForBbox()

    The first returns a 13x13 image of the PRF evalulated at the requested column and row.
    The second trips the 13x13 PRF to match the input bounding box. This facilitiates matching
    the PRF to a Target Pixel File (TPF)


    Notes
    --------
    * API requires the sector to be input as well as the ccd and camera. At the time of writing
      the same model is applicable to all sectors. If this ever changes  the function `sectorLookup()`
      will need to be changed.
    * TODO: For speed, interpolate prfObjects before extracting out 13x13 regular arrays

    """

    def __init__(self, path):
        AbstractPrfLookup.__init__(self, path)
        self.gridSize = 9


    def getPrfForBbox(self, col, row, ccd, camera, sector, bboxIn):
        """Get PRF for a bounding box.

        See `getPrfAtColRow()` and documentation in the same method in the parent class
        """
        args = [ccd, camera, sector]

        return self.abstractGetPrfForBbox(col, row, bboxIn, self.getPrfAtColRow, *args)


    def getPrfAtColRow(self, col, row, ccd, camera, sector):
        """Lookup a 13x13 PRF image for a single location

        Inputs
        ---------
        col, row
            (floats) Location on CCD to lookup. The origin of the CCD is the bottom left.
            Increasing column increases the "x-direction", and row increases the "y-direction"
        ccd
            (int) CCD number. There are 4 CCDs per camera
        camera
            (int) Camera number. The instrument has 4 cameras
        sector
            (int) Sector of observaton.


        Returns
        ---------
        A 13x13 numpy image array.
        """
        col = float(col)
        row = float(row)

        self.checkOutOfBounds(col, row)
        #Currently, the same PRF model applies to all sectors, so
        #we pin the sector number. If the PRF is recalculated at a later
        #date we'll need some logic here.
        sector = self.sectorLookup(sector)
        key = "%1i-%1i-%02i" %(ccd, camera, sector)

        if key not in self.cache:
            self.cache[key] = self.readPrfFile(ccd, camera, sector)

        prfObj = self.cache[key]
        prfArray, evalCols, evalRows = self.getRegularlySampledBracketingPrfs(prfObj, col,row)
        bestPrf = self.interpolatePrf(prfArray, \
            col, row, evalCols, evalRows)

#        regPrf = self.getSingleRegularlySampledPrf(bestPrf, col, row)
        return bestPrf


    def checkOutOfBounds(self, col, row):
        if col < 45 or col > 2091:
            raise ValueError("Requested column (%i) not on phyiscal CCD [45,2091]" %(col))

        if row < 1 or row > 2047:
            raise ValueError("Requested row (%i) not on phyiscal CCD [0,2047]" %(row))


    def sectorLookup(self, sector):
        """Map sector of observation to PRF sector file number.

        At the start of the mission, the same PRFs apply to all sectors. In the future,
        a second PRF file may be released. If that happens, encode the logic of mapping
        sector to file number in this method.

        And remember to update the docstring when you do.

        """
        return 1


    def getRegularlySampledBracketingPrfs(self, prfObj, col, row):
        """Find the 4 grid locations in the PRF file that bracket the requested col,row

        This is an internal function to the class, not intended to be called directly.

        Inputs
        -----------
        prfObj
            (np array). See `readPrfFile()`
        col, row
            (floats) Column and row of interest

        Returns
        ------------
        regPrfArr
            An array of 4 regularly sampled PRFs (regularly sampled PRFs are 13x13 images that can be
            directly compared to a real image (unlike the PRF objects stored on disk, which need to be
            unpacked before use)
        c0
            (np array) Column locations for the 4 images in regPrfArray
        r0
            (np array) Row locations for the 4 images in regPrfArray

        """

        cr = np.array((col, row))


        #Read cols and rows at which PRF is evaluated
        evalCol = npmap(lambda x: x.ccdColumn, prfObj)
        evalRow = npmap(lambda x: x.ccdRow, prfObj)

        #Sort prfArray so its the right shape for getBracketingIndices()
        nEval = np.sqrt(len(prfObj))
        assert nEval == int(nEval), "PRF grid is not square"
        nEval = int(nEval)

        evalColrow = np.vstack((evalCol, evalRow)).transpose()
        srt = np.lexsort((evalRow, evalCol))
        assert isinstance(srt, np.ndarray), "Possibly Py3 problem"

        evalColRow = evalColrow[srt].reshape((nEval, nEval, 2))
        prfArray = prfObj[srt].reshape((nEval, nEval))

        whBracket = getBracketingIndices(evalColRow, cr)

        c0, r0 = [], []
        regPrfArr = []
        for itr in range(4):
            i, j = whBracket[itr]
            #Store column and row
            c0.append( evalColRow[ i, j, 0 ] )
            r0.append( evalColRow[ i, j, 1 ] )

            #Check I did all the book keeping correctly
            assert c0[itr] == prfArray[i,j].ccdColumn
            assert r0[itr] == prfArray[i,j].ccdRow

            #Pull out the 13x13 image for this location
            regPrf = self.getSingleRegularlySampledPrf(prfArray[i,j], col, row)
            regPrfArr.append(regPrf)

        #More checks: check the order of the locations is correct for self.interpolatePrf()
        assert c0[0] == c0[2]
        assert c0[1] == c0[3]
        assert r0[0] == r0[1]
        assert r0[2] == r0[3]
        return np.array(regPrfArr), np.array(c0), np.array(r0)


    def getSingleRegularlySampledPrf(self, singlePrfObj, col, row):
        """
        Look up a regularly sampled PRF. Regularly sampled means sampled at the same
        pixel spacing as the real data.

        Inputs
        ----------
        singlePrfObj
            A prf Obj as returned by `readPrfFile()`
        col, row
            (floats) Column and row of interest

        Returns
        ---------
        A 13x13 image as a numpy 2d array

        Todo
        --------
        This function returns the PRF at the closest point of evaluation. It
        really should interpolate between points to get a PRF that varies
        more smoothly with intrapixel location.
        """

        colOffset, rowOffset = self.getOffsetsFromPixelFractions(singlePrfObj, col, row)
        img = self.getRegularlySampledPrfByOffset(singlePrfObj, colOffset, rowOffset)
        return img


    def getOffsetsFromPixelFractions(self, singlePrfObj, col, row):
        """Private function of `getSingleRegularlySampledPrf()`

        Map the fractional part of the col,row position to an offset into the
        full prf image. For example, if (col, row) = (123,4, 987.6), then
        (colFrac, rowFrac) = (.4, .6).

        This function was developed through trial and error, rather than by
        reference to any design document.
        """
        gridSize = self.gridSize

        colFrac = np.remainder(float(col), 1)
        rowFrac = np.remainder(float(row), 1)

        colOffset = gridSize - np.round(gridSize * colFrac) - 1
        rowOffset = gridSize - np.round(gridSize * rowFrac) - 1

        return int(colOffset), int(rowOffset)


    def getRegularlySampledPrfByOffset(self, singlePrfObj, colOffset, rowOffset):
        """Private function of `getSingleRegularlySampledPrf()`

        The 13x13 pixel PRFs on at each grid location are sampled at a 9x9 intra-pixel grid, to
        describe how the PRF changes as the star moves by a fraction of a pixel in row or column.
        To extract out a single PRF, you need to address the 117x117 array in a slightly funny way
        (117 = 13x9),

        .. code-block:: python

            img = array[ [colOffset, colOffset+9, colOffset+18, ...],
                         [rowOffset, rowOffset+9, ...] ]

        """
        #TODO Reference documentation about samplesPerPixel

        gridSize = self.gridSize

        #TODO: Bounds checks like this don't belong in an inner loop. Remove them
        #once you're sure they never get triggered.
        if colOffset >= gridSize:
            raise ValueError("Requested column offset (%i) too large" %(colOffset))
        if rowOffset >= gridSize:
            raise ValueError("Requested row offset (%i) too large" %(colOffset))


        assert colOffset < gridSize
        assert rowOffset < gridSize

        fullImage = singlePrfObj.values/ float(singlePrfObj.samplesPerPixel)

        #Number of pixels in regularly sampled PRF. Typically 13x13
        nColOut, nRowOut = fullImage.shape
        nColOut /= float(gridSize)
        nRowOut /= float(gridSize)

        iCol = colOffset + (np.arange(nColOut) * gridSize).astype(np.int)
        iRow = rowOffset + (np.arange(nRowOut) * gridSize).astype(np.int)

        #Don't understand why this must be a twoliner
        tmp = fullImage[iRow, :]
        return tmp[:,iCol]


    def interpolatePrf(self, regPrfArray, col, row, evalCols, evalRows):
        """Interpolate between 4 images to find the best PRF at col, row

        This is a private function of the class.

        TODO: Make sure this is right.
        """
        #OK, this is paranoia. These statements have already been asserted, but in a
        #different function
        assert evalCols[0] == evalCols[2]
        assert evalCols[1] == evalCols[3]
        assert evalRows[0] == evalRows[1]
        assert evalRows[2] == evalRows[3]

        p11, p21, p12, p22 = regPrfArray
        c0, c1 = evalCols[:2]
        r0, r1 = evalRows[1:3]

        assert c0 != c1
        assert r0 != r1

        dCol = (col-c0) / (c1-c0)
        dRow = (row-r0) / (r1 - r0)

        #Intpolate across the rows
        tmp1 = p11 + (p21 - p11) * dCol
        tmp2 = p12 + (p22 - p12) * dCol

        #Interpolate across the columns
        out = tmp1 + (tmp2-tmp1) * dRow
        return out


    def readPrfFile(self, ccd, camera, sector):

        if sector > 5:
            raise ValueError("Code needs to be adapted for sectors 6 and above")

        fn = "tess*-00072_035-%i-%i-characterized-prf.mat" %(camera, ccd)
        path = os.path.join(self.path, fn)
        path = glob(path)[0]

        obj = spio.matlab.loadmat(path, struct_as_record=False, squeeze_me=True)
        prfObj = obj['prfStruct']
        return prfObj



def getBracketingIndices(evalColRow, cr):
    """
    Get the indices of `evalColRow` that bracket `cr`

    This is a special function used by TessPrf

    This function encapsulates some fairly knotty bookkeeping. Unless something
    is broken you probably want to leave this function well alone

    Inputs
    --------
    evalColRow
        (3d np array) See discussion below
    cr
        (2 element np array) The column and row to be bracketed

    Returns
    ----------
    A 4x2 numpy array. Each row represents the indices into
    `evalColRow[,,:]` representing the 4 points in `evalColRow`
    that bracket the location represented by evalColRow


    Note
    -----
    The model prf is evaluated on a regular grid across the CCD. Each
    grid point can be represented in two coordinate systems; the
    CCD pixel coordinates (this PRF is evaluated at col,row=24,36,
    and a grid Index (this is the second grid location in column, and
    the third in row). `evalColRow` encodes a mapping from one coord sys
    to the other.

    The zeroth dimension of `evalColRow` encodes the column of the grid
    location (e.g. 2 in the example above). The first dimension
    encodes row of the grid location (3 in the example), the second
    dimension encodes whether the value represents CCD column
    (`evalColRow[:,:,0]`) or CCD row (`evalColRow[:,:,1]`). The
    value at each array element represents the CCD position (either
    column or row).

    The return value of this function is a list of the 4 grid locations
    that bracket the input `cr` in column and row (below left, below right,
    above left, above right)

    Example
    ---------
    `evalColRow` consists of 4 points at which the model prf is evalulated

    .. code-block:: python

        a[0,0,0] =  45
        a[0,0,1] =   1   #Zeroth prf evalulated at (col, row) = (45,1)
        a[0,1,0] =  45
        a[0,1,1] = 128

        a[1,0,0] = 183
        a[1,0,1] =   1
        a[1,1,0] = 183
        a[1,1,1] = 128

        cr = (45, 45)  #Somewhere in the middle

    The return value is

    .. code-block:: python

        [ [0,0], [1,0], [1,0], [1,1] ]

    Because these are the indices that bracket the input col,row
    """
    tmp = (evalColRow - cr)
    dist = np.hypot(tmp[:,:,0], tmp[:,:,1])
    wh = np.unravel_index( np.argmin(dist), dist.shape)

    nearestEval = evalColRow[wh]
    delta = cr - nearestEval

    #Find the 3 other evaluations of the PRF that bracket (col, row)
    tmp = []
    if delta[0] >= 0 and delta[1] >= 0:        #A
        tmp.append( wh )
        tmp.append( wh + np.array((+1, +0)) )
        tmp.append( wh + np.array((+0, +1)) )
        tmp.append( wh + np.array((+1, +1)) )

    elif delta[0] < 0 and delta[1] >= 0:       #S
        tmp.append( wh + np.array((-1, +0)) )
        tmp.append( wh )
        tmp.append( wh + np.array((-1, +1)) )
        tmp.append( wh + np.array((+0, +1)) )

    elif delta[0] < 0 and delta[1] < 0:        #T
        tmp.append( wh + np.array((-1, -1)) )
        tmp.append( wh + np.array((-0, -1)) )
        tmp.append( wh + np.array((-1, +0)) )
        tmp.append( wh )

    else:                                      #C
        tmp.append( wh + np.array((-0, -1)) )
        tmp.append( wh + np.array((+1, -1)) )
        tmp.append( wh )
        tmp.append( wh + np.array((+1, +0)) )
    tmp = np.array(tmp)

    #Check the order of values is correct
    c0 = tmp[:,0]
    r0 = tmp[:,1]
    assert c0[0] == c0[2]
    assert c0[1] == c0[3]
    assert r0[0] == r0[1]
    assert r0[2] == r0[3]

    #Bounds checking
    assert np.min(tmp) >= 0
    assert np.max(tmp[:,0]) < evalColRow.shape[0]
    assert np.max(tmp[:,1]) < evalColRow.shape[1]

    return tmp
