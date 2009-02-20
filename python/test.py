import unittest
import sys, glob, os

sys.path.append( os.path.abspath("bin/gcc-4.1.2/debug" ))

# Do NOT import alignlib here, otherwise there will be weird 
# interactions between tests.
# If alignlib is imported here, it seems than that all 
# tests use the same alignlib.
# Also it seems that tests are all created before being run.
# Thus setting default objects in the constructor will
# affect all subsequent tests.

def suite():
    modules_to_test = map( lambda x: x[:-3], glob.glob("tests/*.py" ))
    
    alltests = unittest.TestSuite()
    for module in map(__import__, modules_to_test):
        alltests.addTest(unittest.findTestCases(module))
    return alltests

if __name__ == '__main__':
    unittest.main(defaultTest='suite')