
import numpy as np

__version__ = "$Id: clipboard.py 2130 2015-09-11 16:56:55Z fmullall $"
__URL__ = "$URL: svn+ssh://flux/home/fmullall/svn/kepler/k2phot/clipboard.py $"


class Clipboard(object):
    """A clipboard for holding a series of key value pairs

    A value can be any valid python object. A clipboard
    is very similar to a normal python dictionary except

    * Sub dictionaries are easier to address.
      x = Clipboard()
      x['y'] = dict('a':1)
      print x['y.a']  #Returns x['y']['a]

    * When accessing numpy arrays, a copy is returned by
      default, which is usually what I want.

    * Printing the contents of the clip is much prettier than
      for a dictionary.
    """

    def __init__(self, *args, **kwargs):
        """
        Input:
        --------
        A dictionary, or any allowed initialisation arguments for
        a dictionary
        """
        self.store = dict(*args, **kwargs)

    def __getitem__(self, keyString):
        return self.get(keyString)

    def __setitem__(self, keyString, value):
        self.set(keyString, value)

    def __repr__(self):
        return self.asString()

    def __delitem__(self, keyString):
        keys = keyString.split(".")

        tmp = self.store
        for k in keys[:-1]:
            tmp = tmp[k]

        del tmp[keys[-1]]


#    def __dir__(self):
#        return self.getFullKeyList()


    def getFullKeyList(self, arg=None, prefix=None):
        if arg is None:
            arg = self.store
        if prefix is None:
            prefix = ""
        else:
            prefix = "%s." %(prefix)

        keyList = []

#        import pdb; pdb.set_trace()
        for k in arg:
            value = arg[k]
            if isinstance(value, dict) or isinstance(value, Clipboard):
                prefix = "%s%s" %(prefix, k)
                keyList.extend(self.getFullKeyList(value, prefix))
            else:
                keyList.append("%s%s" %(prefix, k))
#            print k, keyList
        return keyList


    def __getattr__(self, keyString):
        return self.get(keyString)


    def unsetException(self):
        if 'exception' in self.store:
            del self.store['exception']

    def get(self, keyString, defaultValue=None, npCopy=True):
        """Access a value.

        Inputs:
        ---------
        keyString
            (string) A dictionary key, or a set of keys joined by
            full stops.

        Optional Inputs:
        ------------------
        defaultValue
            (default **None**). If requested key does not exist, return
            this value instead. If the defaultValue is **None** raise
            a KeyError instead

        npCopy
            (default: **True**) Return a copy of numpy arrays instead
            of a pointer to one. This prevents the most common corruption
            of Clipboard contents.


        Returns:
        ---------
        Contents of the clipboard()

        Example:
        ----------
        ::

            clip = Clipboard()
            clip['a']['b']  #Returns the same as
            clip['a.b']

        you can also write
        ::

            clip.get('a.b', defaultValue=0)

        """
        keys = keyString.split(".")

        tmp = self.store
        for k in keys:
            try:
                tmp = tmp[k]
            except KeyError:
                if defaultValue is not None:
                    return defaultValue

                raise KeyError("Required key %s not found while looking up %s" %(k, keyString))

        if isinstance(tmp, np.ndarray):
            if npCopy:
                return tmp.copy()
        return tmp


    def set(self, keyString, value):
        """Set a value in the Clipboard.

        See doc string for get() for details

        Returns:
        ----------
        **None**
        """
        keys = keyString.split(".")

        tmp = self.store
        for k in keys[:-1]:
            try:
                tmp = tmp[k]
            except KeyError:
                raise KeyError("Required key %s not found while looking for %s" %(k, keyString))
        tmp[keys[-1]] = value


    def keys(self):
        return self.store.keys()


    def getSubclip(self, keyString, defaultValue=None, npCopy=True):
        """Same as get(), but if the returned value is a dictionary
        upcast it to a Clipboard()
        """

        val = self.get(keyString, defaultValue, npCopy)
        if isinstance(val, dict):
            return Clipboard(val)


    def pprint(self, prefix="", maxLevel=np.inf):
        """Pretty print the contents of a clip

        Optional Inputs:
        ----------------
        prefix
            Don't change this value, it's used to track the
            level of recursion

        Returns:
        ---------
        **None**

        Output:
        ---------
        Prints text to stdout.
        """
        print self.asString(prefix, maxLevel)


    def asString(self, prefix="", maxLevel=np.inf):
        """Pretty print the contents of a clip

        Optional Inputs:
        ----------------
        prefix
            Don't change this value, it's used to track the
            level of recursion

        maxLevel
            Only descend this many levels into the clipboard
        Returns:
        ---------
        A string

        """

        out = []
        clip = self.store
        for k in sorted(clip.keys()):
            label = "%s%s:" %(prefix, k)
            v = clip[k]
            if isinstance(v, dict):
                v = Clipboard(v)

            if isinstance(v, Clipboard):
                out.append(label)
                if len(prefix) < maxLevel-1:
                    out.append(v.asString(prefix+"+", maxLevel))
                else:
                    out[-1] = "%s Clip" %(out[-1])
            elif isinstance(v, np.ndarray):
                out.append("%s np.ndarray %s" %(label, str(v.shape)))
            elif isinstance(v, str):
                mx = min(50, len(v))
                out.append("%s %s" %(label, v[:mx]))
            elif isinstance(v, int):
                out.append("%s %i" %(label, v))
            elif isinstance(v, float):
                out.append("%s %g" %(label, v))
            elif isinstance(v, tuple):
                out.append("%s %g" %(label, str(v)))
            elif isinstance(v, list):
                mx = min(10, len(v))
                out.append("%s %s" %(label, v[:mx]))

        return "\n".join(out)



def test():
    c = Clipboard()
    c['x'] = dict()
    c['x.y'] = dict()
    c['x.y.z'] = "A str"

    print c.asString(maxLevel=2)

    del c['x.y.z']

    c['exception'] = True
    c.unsetException()
    return c
