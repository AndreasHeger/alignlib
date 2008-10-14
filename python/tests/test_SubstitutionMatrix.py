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

        self.mMatrices = [ (lambda : makeSubstitutionMatrix( 23, 1, -1 ), 23),
                          (makeSubstitutionMatrixBlosum62, 23),
                          (makeSubstitutionMatrixBlosum50, 23),
                          (makeSubstitutionMatrixPam250, 23),
                          (makeSubstitutionMatrixPam120, 23),
                          (lambda : makeSubstitutionMatrixBackTranslation( 1, -1, 0.5, getEncoder( Protein23) ), 128), 
                   ]
                    
    def testMake(self):
        """check if all matrices can be created and are square."""
        for matrix, size in self.mMatrices:
            m = matrix()
            self.assertEqual( m.getNumRows(), m.getNumCols() )
            self.assertEqual( m.getNumRows(), size )    
        
    def testSetDefault(self):
        
        for matrix, size in self.mMatrices:
            setDefaultSubstitutionMatrix( matrix() )
            
    def testGetDefault(self):
        
        matrix = getDefaultSubstitutionMatrix()
        
        
def suite():
    suite = unittest.TestSuite()
    suite.addTest(SubstitutionMatrixTestCase)
    return suite

if __name__ == "__main__":
    unittest.main()


        





