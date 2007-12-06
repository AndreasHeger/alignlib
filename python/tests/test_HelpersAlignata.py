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
import unittest, tempfile
import StringIO, cStringIO
from alignlib import *

class HelpersAlignataTestCase( unittest.TestCase ):

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
                    
    def testWritePairAlignment(self ):
        
        row = makeSequence( "AAAAAAAAAAAAAAAAAAAAAAAAAAAAA" )
        col = makeSequence( "AAAAAAAAAAAAAAAAAAAAAAAAAAAAA" )

        self.buildAlignment( self.mAlignata )
        
        outfile = open("test.out", "wb" )
        writePairAlignment(outfile, row, col, self.mAlignata )
        outfile.close()
        
        outfile = tempfile.TemporaryFile( "w+b" )
        result = writePairAlignment( outfile, row, col, self.mAlignata )
        outfile.seek(0)
        input = outfile.readlines()
            
    def testAlignataCompressed( self ):

        self.buildAlignment( self.mAlignata )

        outfile = tempfile.TemporaryFile( "w+b" )
        writeAlignataCompressed( outfile, self.mAlignata )
        outfile.seek(0)
        row_ali, col_ali = outfile.readline().split("\t")
        
        alignment = self.mAlignata.getNew()
        
        fillAlignataCompressed( alignment,
                                self.mAlignata.getRowFrom(),
                                row_ali, 
                                self.mAlignata.getColFrom(),
                                col_ali )
        
        self.assertTrue( self.compareResidueWise( self.mAlignata, alignment ) )
    
    def compareResidueWise(self, a, b, inverse = False):
        """compare two alignments. Returns true if they are identical
        residuewise."""
        
        
        for x in range( a.getRowFrom(), a.getRowTo() ):
            if a.mapRowToCol( x ) != b.mapRowToCol( x ):
                return False
            
        return True
        
        it1 = a.begin()
        it1_end = a.end()
        
        it2 = b.begin()
        it2_end = b.end()
        
        is_identical = True
        while it1 != it1_end:
            row = it1.getReference()
            col = it2.getReference()
            if not inverse: 
                if it1.mRow != it2.mRow and it1.mCol != it2.mCol:
                    is_identical = False 
            else:
                if it1.mRow != it2.mCol and it1.mCol != it2.mRow:
                    is_identical = False
                    
        return is_identical
        
def suite():
    suite = unittest.TestSuite()
    suite.addTest(AlignataVectorTestCase)
    suite.addTest(AlignataSetTestCase)    
    return suite

if __name__ == "__main__":
    unittest.main()
