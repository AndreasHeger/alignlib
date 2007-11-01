# alignlib - a library for aligning protein sequences
# 
# $Id: test_Alignata.py,v 1.3 2004/01/23 17:34:58 aheger Exp $
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

class AlignataTestCase( unittest.TestCase ):

    def setUp( self ):

        self.mAlignata = makeAlignataVector()

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
        self.buildAlignment( self.mAlignata )
        
    def testGetRowFrom( self ):
        self.testAddPair()
        self.assertEqual( self.mAlignata.getRowFrom(), 3) 

    def testGetRowTo( self ):
        self.testAddPair()
        self.assertEqual( self.mAlignata.getRowTo(), 12) 

    def testGetColFrom( self ):
        self.testAddPair()
        self.assertEqual( self.mAlignata.getRowFrom(), 3) 

    def testGetColTo( self ):
        self.testAddPair()
        self.assertEqual( self.mAlignata.getColTo(), 12) 

    def testRemoveRowRegion( self ):
        self.testAddPair()

        f = self.mAlignata.removeRowRegion
        f( 3, 5)
        self.assertEqual( self.mAlignata.getRowFrom(), 6)
        self.assertEqual( self.mAlignata.getColFrom(), 7)        
        f( 3, 5 )
        f( 10, 12 )
        self.assertEqual( self.mAlignata.getRowTo(), 9)
        self.assertEqual( self.mAlignata.getColTo(), 10)        
        f( 1, 20 )
        self.assertEqual( self.mAlignata.getLength(), 0 )
        self.assertEqual( self.mAlignata.isEmpty(), True)

        ## try nonsense input
        self.testAddPair()
        f(-1, -3)
        f(20, 30)
        f(0, 2)
        f(1, 2)
        f(2, 1)                                
        self.assertEqual( self.mAlignata.getLength(), 12 )

    def testRemoveColRegion( self ):
        self.testAddPair()

        f = self.mAlignata.removeColRegion
        f( 3, 5)
        self.assertEqual( self.mAlignata.getRowFrom(), 5)
        self.assertEqual( self.mAlignata.getColFrom(), 6)        
        f( 3, 5 )
        f( 10, 12 )
        self.assertEqual( self.mAlignata.getRowTo(), 8)
        self.assertEqual( self.mAlignata.getColTo(), 8)        
        f( 1, 20 )
        self.assertEqual( self.mAlignata.getLength(), 0 )
        self.assertEqual( self.mAlignata.isEmpty(), True)

        ## try nonsense input
        self.testAddPair()
        f(-1, -3)
        f(20, 30)
        f(0, 2)
        f(1, 2)
        f(2, 1)                                
        self.assertEqual( self.mAlignata.getLength(), 12 )

    def testSwitchRowCol( self ):
        """test switching of row and column and mapping."""
        
        self.testAddPair()
        ali = self.mAlignata.getClone()
        ali.switchRowCol()

        for x in range(self.mAlignata.getRowFrom(), self.mAlignata.getRowTo() + 1):
            y = self.mAlignata.mapRowToCol( x )
            if y:
                z = ali.mapRowToCol( y )
                self.assertEqual( x, z )
            y = self.mAlignata.mapColToRow( x )
            if y:
                z = ali.mapColToRow( y )
                self.assertEqual( x, z )

    def testIterator( self ):
        '''!!!TODO!!!'''        
        self.testAddPair()
                
class AlignataVectorTestCase( AlignataTestCase ):

    def setUp( self ):
        self.mAlignata = makeAlignataVector()
        
class AlignataSetTestCase( AlignataTestCase ):

    def setUp( self ):
        self.mAlignata = makeAlignataSet()

class AlignataSetColTestCase( AlignataTestCase ):

    def setUp( self ):
        self.mAlignata = makeAlignataSetCol()

class AlignataHashTestCase( AlignataTestCase ):

    def setUp( self ):
        self.mAlignata = makeAlignataHash()


class AlignataHashDiagonalTestCase( AlignataTestCase ):
    '''!!!TODO!!! fix me'''
    def setUp( self ):
        # self.mAlignata = makeAlignataHashDiagonal()
        self.mAlignata = makeAlignataVector()

class AlignataMatrixRowTestCase( AlignataTestCase ):
    '''!!!TODO!!! fix me'''
    def setUp( self ):
        # self.mAlignata = makeAlignataMatrixRowDiagonal()
        self.mAlignata = makeAlignataVector()

class AlignataMatrixUnsortedTestCase( AlignataTestCase ):
    '''!!!TODO!!! fix me'''
    def setUp( self ):
        # self.mAlignata = makeAlignataMatrixUnsorted()
        self.mAlignata = makeAlignataVector()

class AlignataMatrixDiagonalTestCase( AlignataTestCase ):
    '''!!!TODO!!! fix me'''
    def setUp( self ):
        # self.mAlignata = makeAlignataMatrixDiagonal()
        self.mAlignata = makeAlignataVector()

def suite():
    suite = unittest.TestSuite()
    suite.addTest(AlignataVectorTestCase)
    suite.addTest(AlignataSetTestCase)    
    return suite

if __name__ == "__main__":
    unittest.main()
