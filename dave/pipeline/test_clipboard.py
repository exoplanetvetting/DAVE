# -*- coding: utf-8 -*-
"""
Created on Fri Dec  4 14:24:06 2015

@author: fergal

$Id$
$URL$
"""

__version__ = "$Id$"
__URL__ = "$URL$"



from clipboard import Clipboard, loadClipboard
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


    def testGetAttributes(self):
        """Doesn't work yet"""
        c = self.clip
        self.assertEqual(c.c2.a, 1)
        self.assertEqual(c.c3.y, 2)


    #def testSetAttributes(self):
        #"""This test fails. Group is assigned an attribute of the
        #class, instead of being cast as a clipboard and stored in
        #object.store['group']. This gives very confusing results"""
        #c = self.clip
        #c.array = np.zeros((10))
        #c.group = dict()
        #c.group.x = 1

        #self.assertEqual(len(c.array), 10)
        #self.assertEqual(c.group.x, 1)

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


    def test_get(self):
        x = self.clip.get('c2.NotAKey', 5)
        self.assertEquals(x, 5)

#Saveing from clipboard turned off right now
#    def test_save(self):
#        self.clip.save('tmp.shelf')
#        c2 = loadClipboard("tmp.shelf")
#
#        self.assertEqual(c2['c2.a'], 1)
#        self.assertEqual(c2['c3.x'], 1)

    def testInit(self):
        c = Clipboard(a=1, b=2, c={'x':1, 'y':2})
        self.assertEqual(c['a'], 1)
        self.assertEqual(c['b'], 2)
        self.assertEqual(c['c.x'], 1)
        self.assertEqual(c['c.y'], 2)

        d={'x':1, 'y':2}
        c = Clipboard(a=d)
        self.assertEqual(c['a.x'], 1)
        self.assertEqual(c['a.y'], 2)

        c = Clipboard(d)
        self.assertEqual(c['x'], 1)
        self.assertEqual(c['y'], 2)

    def testAddedDictionaryBecomesClipboard(self):
        c = self.clip
        d = {'x':1, 'y':2}

        c['d'] = d
        self.assertTrue(isinstance(c.d, Clipboard))
#        self.assert(True)
        self.assertEqual(c.d.x, 1)
        self.assertEqual(c.d.y, 2)

    def testAddingTuple(self):
        """Reproducing Issue #7. Adding a tuple made asString() fail"""
        c = self.clip

        c['aTuple'] = ('a', 'b', 'c')
        c.asString()

if __name__ == "__main__":
    unittest.main()

