import unittest
import sys, glob, os

sys.path.append( os.path.abspath("bin/gcc-4.1.2/debug" ))

import alignlib

def suite():
    modules_to_test = map( lambda x: x[:-2], glob.glob("tests/*.py" )) # and so on
    alltests = unittest.TestSuite()
    for module in map(__import__, modules_to_test):
        alltests.addTest(unittest.findTestCases(module))
    return alltests

if __name__ == '__main__':
    unittest.main(defaultTest='suite')