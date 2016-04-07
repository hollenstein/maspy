from future.utils import viewkeys, viewvalues, viewitems, listvalues, listitems
import unittest

import maspy.new
import numpy
import tempfile

class Test_si(unittest.TestCase):
    def test_1(self):
        array1 = numpy.array(range(1000), dtype=numpy.float64)
        array2 = numpy.array(range(100, 2000, 2), dtype=numpy.float64)
        adict = {'1': array1, '2': array2}

        with tempfile.NamedTemporaryFile() as filedummy:
            maspy.new.saveToFile(adict, filedummy)
            filedummy.seek(0)

            bdict = maspy.new.loadFromFile(filedummy)
            self.assertSetEqual(set(viewkeys(adict)), set(viewkeys(bdict)))
            for key in adict:
                self.assertTrue(numpy.array_equal(adict[key], bdict[key]))
