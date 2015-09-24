# -*- coding: utf-8 -*-
import numpy as np

__version__ = "$Id$"
__URL__ = "$URL$"

class Nca(np.ndarray):
    """
        A thin wrapper around numpy to allow column addressing by name
        instead of number.

        I frequently have a large 2d array where each row represents an
        object, and each column represent an attribute of that object.
        For example::
          Obj1   Period  Epoch Radius
          Obj2   Period  Epoch Radius

        To get all the radii, I'd like to say
        array[:, "radius"] instead of array[:,3]
        This class allows me to do that. Here is an example::
          row = '0 1 2 3'.split()
          col = 'a b c d'.split()
          nameDict = dict()
          nameDict[0] = row
          nameDict[1] = col

          data = np.arange(16).reshape(4,4) + 1
          ca = Nca( data , nameDict)

        Anywhere you normally put a number, you can now use a string.
        For example::
          ca['0', 'b']
          ca['0', :'b']
          ca[idx, 'a':'c']  etc.

        The return behaviour is similar to numpy, except that an object
        of the extended class is returned where numpy would return an array.
        The lookup information is also transferred in a natural manner.
        For example::
          x = ca[:, 'c']
          x[:, 'b'] #Doesn't work
          x[:, 'c'] #returns x[:,0]


        Inputs:
        ---------
        ndArray
            A numpy array
        nameDict
            A lookup table matching strings to array positions.
            See setLookup() for more details


        Notes:
        ------
        * How to subclass numpy's ndarray is taken from
          <http://docs.scipy.org/doc/numpy/user/basics.subclassing.html>_
          The __new__ and __array_finalize__ methods are copied from there.

        * There is one known bug with this code.arr[:4] doesn't return
          the correct lookup table. As a workaround
          use arr[:,4, :], which does.

        Todo:
        -------
        * Fix bug where array[:4] fails (but array[:'4'] doesn't
        * Add a metadata dictionary
        * Write tests to make sure lookup dictionary is being correctly sliced

    """


    def __new__(cls, input_array, nameDict=None):
        obj = np.asarray(input_array).view(cls)
        obj.lookup = nameDict
        return obj

    def __array_finalize__(self, obj):
        if obj is None:
            return
        self.lookup = getattr(obj, 'lookup', None)

    def __getitem__(self, key):
        #print "Input key", key, type(key)
        key = self.parseKey(key)
        #print "Parsed key", key
        returnObj = np.ndarray.__getitem__(self, key)

        #
        #This stuff fails for arr[:4]
        #
        if isinstance( returnObj, Nca):
            if self.lookup is None:
                newLookup = None
            else:
                newLookup = self._setNewLookup(key)
            return Nca(returnObj, nameDict=newLookup)
        else:
            return returnObj


    def __setitem__(self, key, value):
        key = self.parseKey(key)
        np.ndarray.__setitem__(self, key, value)


    def _setNewLookup(self, key):

        #If only 1 int, we're looking at one row of zeroth dimension.
        if isinstance(key, int):
            try:
                return self.lookup[0][key]
            except KeyError:
                return None

        #Similarly for a slice
        if isinstance(key, slice):
            return self.lookup[0][key]

        #If it's a tuple, we have multiple dimensions
        if isinstance(key, tuple):
            newLookup = dict()
            for i in range(len(key)):
                try:
                    newLookup[i] = np.array(self.lookup[i])[ key[i]]
                except KeyError:
                    #Lookup not defined for this dimension
                    continue

                newLookup[i] = list(newLookup[i])
            return newLookup


    def parseKey(self, key, dim=0):
        #import pdb; pdb.set_trace()
        if isinstance(key, str):
            try:
                key = self.lookup[dim].index(key)
            except ValueError:
                raise KeyError(\
                    "key '%s' not a recognised column in dimension %i" %(key, dim))

        if isinstance(key, list):
            #raise NotImplementedError("List should be supported but aren't")
            tmp = key
            for i in  range(len(tmp)):
                tmp[i] = self.parseKey(tmp[i], dim=dim)
            #import pdb; pdb.set_trace()
            return tmp

        if isinstance(key, tuple):
            tmp = list(key)
            for i in range(len(tmp)):
                tmp[i] = self.parseKey(tmp[i], dim=i)
            return tuple(tmp)

        if isinstance(key, slice):
            start = self.parseKey(key.start, dim=dim)
            stop = self.parseKey(key.stop, dim=dim)
            step = key.step
            return slice(start, stop, step)

            return key

        #No more strings to strip out
        return key


    def setLookup(self, dim, colNameList):
        """Set the mapping from column name to column number.

        Inputs:
        dim             (int) Which dimension to create a lookup table for.
        colNameList     (list) What to call the column names. colNameList
                        must be the same length as self.shape[dim]


        For example
        setLookup(1, ['time', 'flux', 'unc'])

        allows you to access an array with
        x[:, 'time']
        """

        if not isinstance(dim, int):
            raise TypeError("Dim must be an integer")

        if len(colNameList) != self.shape[dim]:
            raise ValueError(" Length of colNameList (%i) not equal length of array in dimension %i (%i)" \
                %(len(colNameList), dim, self.shape[dim]))

        if self.lookup is None:
            self.lookup = dict()
        self.lookup[dim] = list(colNameList)


    def asarray(self):
        return self.view(np.ndarray)



def example():
    row = '0 1 2 3'.split()
    col = 'a b c d'.split()
    nameDict = dict()
    nameDict[0] = row
    nameDict[1] = col

    data = np.arange(16).reshape(4,4) + 1
    ca = Nca( data , nameDict)
    return ca

