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
import unittest
from alignlib import *

class AlignmentTestCase( unittest.TestCase ):

    def setUp( self ):

        self.mAlignment = makeAlignmentVector()

    def buildAlignment( self, alignment ):
        
        alignment.clear()
        
        for d in (
            (3,3,1),    
            (4,4,1),
            (5,6,1),
            (6,7,1),
            (8,8,1),
            (9,10,1),
            (10,11,1),
            (12,12,1)):
            alignment.addPair( d[0], d[1], float(d[2]))

    def testAddPair( self ):
        self.buildAlignment( self.mAlignment )
        
    def testGetRowFrom( self ):
        self.testAddPair()
        self.assertEqual( self.mAlignment.getRowFrom(), 3) 

    def testGetRowTo( self ):
        self.testAddPair()
        self.assertEqual( self.mAlignment.getRowTo(), 13) 

    def testGetColFrom( self ):
        self.testAddPair()
        self.assertEqual( self.mAlignment.getRowFrom(), 3) 

    def testGetColTo( self ):
        self.testAddPair()
        self.assertEqual( self.mAlignment.getColTo(), 13) 

    def testRemoveRowRegion( self ):
        self.testAddPair()

        f = self.mAlignment.removeRowRegion
        f( 3, 5)
        self.assertEqual( self.mAlignment.getRowFrom(), 5)
        self.assertEqual( self.mAlignment.getColFrom(), 6)        
        f( 10, 13 )
        self.assertEqual( self.mAlignment.getRowTo(), 10)
        self.assertEqual( self.mAlignment.getColTo(), 11)        
        f( 0, 20 )
        self.assertEqual( self.mAlignment.getLength(), 0 )
        self.assertEqual( self.mAlignment.isEmpty(), True)

        ## try nonsense input
        self.testAddPair()
        f(-1, -3)
        f(20, 30)
        f(0, 2)
        f(1, 2)
        f(2, 1)                                
        self.assertEqual( self.mAlignment.getLength(), 12 )

    def testRemoveColRegion( self ):
        self.testAddPair()

        f = self.mAlignment.removeColRegion
        f( 3, 5)
        self.assertEqual( self.mAlignment.getRowFrom(), 5)
        self.assertEqual( self.mAlignment.getColFrom(), 6)        
        f( 10, 13 )
        self.assertEqual( self.mAlignment.getRowTo(), 9)
        self.assertEqual( self.mAlignment.getColTo(), 9)        
        f( 0, 20 )
        self.assertEqual( self.mAlignment.getLength(), 0 )
        self.assertEqual( self.mAlignment.isEmpty(), True)

        ## try nonsense input
        self.testAddPair()
        f(-1, -3)
        f(20, 30)
        f(0, 2)
        f(1, 2)
        f(2, 1)                                
        self.assertEqual( self.mAlignment.getLength(), 12 )

    def testSwitchRowCol( self ):
        """test switching of row and column and mapping."""
        
        self.testAddPair()
        ali = self.mAlignment.getClone()
        ali.switchRowCol()

        for x in range(self.mAlignment.getRowFrom(), self.mAlignment.getRowTo()):
            y = self.mAlignment.mapRowToCol( x )
            if y >= 0:
                z = ali.mapRowToCol( y )
                self.assertEqual( x, z )
                
            y = self.mAlignment.mapColToRow( x )            
            if y >= 0:
                z = ali.mapColToRow( y )
                self.assertEqual( x, z )

    def testIterator( self ):
        '''!!!TODO!!!'''        
        self.testAddPair()
                
class AlignmentVectorTestCase( AlignmentTestCase ):

    def setUp( self ):
        self.mAlignment = makeAlignmentVector()

class AlignmentSetTestCase( AlignmentTestCase ):

    def setUp( self ):
        self.mAlignment = makeAlignmentSet()

class AlignmentSetColTestCase( AlignmentTestCase ):

    def setUp( self ):
        self.mAlignment = makeAlignmentSetCol()

class AlignmentHashTestCase( AlignmentTestCase ):

    def setUp( self ):
        self.mAlignment = makeAlignmentHash()


class AlignmentHashDiagonalTestCase( AlignmentTestCase ):
    '''!!!TODO!!! fix me'''
    def setUp( self ):
        self.mAlignment = makeAlignmentHashDiagonal()
        self.mAlignment = makeAlignmentVector()

class AlignmentMatrixRowTestCase( AlignmentTestCase ):
    '''!!!TODO!!! fix me'''
    def setUp( self ):
        self.mAlignment = makeAlignmentMatrixRowDiagonal()
        self.mAlignment = makeAlignmentVector()

class AlignmentMatrixUnsortedTestCase( AlignmentTestCase ):
    '''!!!TODO!!! fix me'''
    def setUp( self ):
        self.mAlignment = makeAlignmentMatrixUnsorted()
        self.mAlignment = makeAlignmentVector()

class AlignmentMatrixDiagonalTestCase( AlignmentTestCase ):
    '''!!!TODO!!! fix me'''
    def setUp( self ):
        self.mAlignment = makeAlignmentMatrixDiagonal()
        self.mAlignment = makeAlignmentVector()

def suite():
    suite = unittest.TestSuite()
    suite.addTest(AlignmentVectorTestCase)
    suite.addTest(AlignmentSetTestCase)    
    return suite

if __name__ == "__main__":
    unittest.main()
