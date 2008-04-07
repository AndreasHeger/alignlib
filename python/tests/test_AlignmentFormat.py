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

class AlignmentFormatTestCase( unittest.TestCase ):

    def setUp( self ):
        self.mAlignment = makeAlignmentVector()
        addDiagonal2Alignment( self.mAlignment, 5, 10, 0)
        addDiagonal2Alignment( self.mAlignment, 10, 15, 5)
        addDiagonal2Alignment( self.mAlignment, 25, 30, -5)
        self.mFormat = None

    def testCopy(self):
        
        if self.mFormat:
            a = str(self.mFormat( self.mAlignment ))
            alignment = makeAlignmentVector()
            self.mFormat(a).copy( alignment )
            b = str(self.mFormat( alignment ))
            self.assertEqual( a, b )
            
    def testLoad(self):

        if self.mFormat:
            a = str(self.mFormat( self.mAlignment ))
            alignment = makeAlignmentVector()
            x = self.mFormat()
            x.load( a)
            b = str(x)
            self.assertEqual( a, b )

class AlignmentFormatBlocksTestCase( AlignmentFormatTestCase ):

    def setUp(self):
        AlignmentFormatTestCase.setUp( self )
        self.mFormat = AlignmentFormatBlocks
        
class AlignmentFormatEmissionsTestCase( AlignmentFormatTestCase ):
    
    def setUp(self):
        AlignmentFormatTestCase.setUp( self )
        self.mFormat = AlignmentFormatEmissions

    def testInput2(self):
        a = self.mFormat( 5, "+5-5+10", 5, "+10-5+5" )
        self.assertEqual( a.mRowTo, 20 )
        self.assertEqual( a.mColTo, 20 )
        
class AlignmentFormatDiagonalsTestCase( AlignmentFormatTestCase ):
    
    def setUp(self):
        AlignmentFormatTestCase.setUp( self )
        self.mFormat = AlignmentFormatDiagonals
                
def suite():
    suite = unittest.TestSuite()
    suite.addTest(AlignmentFormatBlocksTestCase)
    suite.addTest(AlignmentFormatEmissionsTestCase)    
    suite.addTest(AlignmentFormatDiagonalsTestCase)    
    return suite

if __name__ == "__main__":
    unittest.main()


        





