# -*- coding: utf-8 -*-
"""
Created on Fri Dec  4 14:24:06 2015

@author: fergal

$Id$
$URL$
"""

__version__ = "$Id$"
__URL__ = "$URL$"



from clipboard import Clipboard
import unittest
import numpy as np


class TestClipboard(unittest.TestCase):

    def setUp(self):
        c2 = Clipboard(a=1, b=2)
        c3 = {'x':1, 'y':2}
        c4 = "A string"

        clip = Clipboard(c2=c2, c3=c3, c4=c4)
        self.clip = clip


    def test_smoke(self):
        self.clip
        self.clip['c2']
        self.clip['c3']
        self.clip['c4']

        self.assertEqual(self.clip['c2.a'], 1)
        self.assertEqual(self.clip['c3.x'], 1)

    def test_Unset(self):
        self.clip['exception'] = True

        self.clip.unsetException()
        with self.assertRaises(KeyError):
            self.clip['exception']

    def test_keyList(self):
        expected = "c2.a c2.b c3.x c3.y c4".split()

        msg = "KeyList: %s" % str(self.clip.getFullKeyList())
        self.assertItemsEqual(self.clip.getFullKeyList(), expected, msg)

    def test_numpyCopy(self):
        c = self.clip
        c['array'] = np.zeros(10)

        tmp = c['array']
        tmp[0] = 10

        self.assertEqual( int(np.sum(c['array'])), 0)

if __name__ == "__main__":
    unittest.main()

