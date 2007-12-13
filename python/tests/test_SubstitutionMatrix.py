# alignlib - a library for aligning protein sequences
# 
# $Id: test_Alignment.py,v 1.3 2004/01/23 17:34:58 aheger Exp $
# 
# Copyright (C) 2004 Andreas Heger
# 
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
#
import unittest, sys
from alignlib import *

class SubstitutionMatrixCase( unittest.TestCase ):

    def setUp( self ):
        
        self.mMatrixFactory = makeSubstitutionMatrixAAIdentity
        
    def testSetDefault(self):
    
        setDefaultSubstitutionMatrix( self.mMatrixFactory() )
        setDefaultSubstitutionMatrix( self.mMatrixFactory() )
        
        
def suite():
    suite = unittest.TestSuite()
    suite.addTest(SubstitutionMatrixTestCase)
    return suite

if __name__ == "__main__":
    unittest.main()


        





