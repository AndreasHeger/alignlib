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
import sys

class AlignatorTestCase( unittest.TestCase ):

    def setUp( self ):

        self.mReferenceSequence1 = "AAAAAAACCCCAAAAAAA"
        self.mReferenceSequence2 = "WWAAAAAWWWWAAAAAWW"         
        self.mAlignataA2B = makeAlignataVector()
        self.mAlignataB2A = makeAlignataVector()        
        self.mSeq1 = makeSequence( self.mReferenceSequence1 )
        self.mSeq2 = makeSequence( self.mReferenceSequence2 )

        ## chose two sequences, so that with -10.0, -1.0 gap penalties
        ## there will be no internal gaps
        self.mPro1 = makeProfile( self.mReferenceSequence1 * 2, 2)
        self.mPro2 = makeProfile( self.mReferenceSequence2 * 2, 2)        

        self.mSeqs = [self.mSeq1, self.mSeq2, self.mPro1, self.mPro2 ]

        self.mSames = [ (0,2), (1,3) ]
            
        self.mPairs = []
        for x in range(len(self.mSeqs)):
            for y in range(x, len(self.mSeqs)):
                self.mPairs.append( (x, y) )

        self.mAlignator = None
        
    def testAlignment( self ):

        if not self.mAlignator: return
        
        for row, col in self.mPairs:
            self.mAlignataA2B.clear()
            self.mAlignataB2A.clear()            
            self.mAlignator.align( self.mSeqs[row], self.mSeqs[col], self.mAlignataA2B )
            self.mAlignator.align( self.mSeqs[col], self.mSeqs[row], self.mAlignataB2A )            
            self.checkAlignment( row, col )

    def checkAlignment( self, row, col ):
        """general sanity checks."""

        self.assertEqual( self.mAlignataA2B.getLength(), self.mAlignataB2A.getLength() )
        self.assertEqual( self.mAlignataA2B.getScore(), self.mAlignataB2A.getScore() )
        self.assertEqual( self.mAlignataA2B.getNumGaps(), self.mAlignataB2A.getNumGaps() )
        self.assertEqual( self.mAlignataA2B.getRowFrom(), self.mAlignataB2A.getColFrom() )
        self.assertEqual( self.mAlignataA2B.getRowTo(), self.mAlignataB2A.getColTo() )        
        
class AlignatorDPGlobalWithEndGapsPenaltiesTestCase( AlignatorTestCase ):

    def setUp( self ):
        AlignatorTestCase.setUp( self )
        # penalize for gaps at the ends. This will force the residue to be aligned
        self.mAlignator = makeAlignatorDPFull( ALIGNMENT_GLOBAL, -10.0, -1.0, True, True )

    def checkAlignment( self, row, col ):

        AlignatorTestCase.checkAlignment( self, row, col )
        self.assertEqual( self.mAlignataA2B.getLength(), self.mSeqs[row].getLength() )            

class AlignatorDPGlobalNoEndGapsPenaltiesTestCase( AlignatorTestCase ):

    def setUp( self ):
        AlignatorTestCase.setUp( self )
        # do not penalize for gaps at the ends (default). This will allow mismatching residues
        # at the ends to be unaligned
        self.mAlignator = makeAlignatorDPFull( ALIGNMENT_GLOBAL, -10.0, -1.0, False, False )

    def checkAlignment( self, row, col ):

        AlignatorTestCase.checkAlignment( self, row, col )
        
        if row == col or (row,col) in self.mSames:
            self.assertEqual( self.mAlignataA2B.getLength(), self.mSeqs[row].getLength() )            
        else:
            self.assertEqual( self.mAlignataA2B.getLength(), self.mSeqs[row].getLength() - 4 )

class AlignatorDPLocalTestCase( AlignatorTestCase ):

    def setUp( self ):
        AlignatorTestCase.setUp( self )
        # penalize for gaps at the ends. This will force the residue to be aligned
        self.mAlignator = makeAlignatorDPFull( ALIGNMENT_LOCAL, -10.0, -1.0 )

    def checkAlignment( self, row, col ):

        AlignatorTestCase.checkAlignment( self, row, col )

        if row == col or (row,col) in self.mSames:
            self.assertEqual( self.mAlignataA2B.getLength(), self.mSeqs[row].getLength() )            
        else:
            self.assertEqual( self.mAlignataA2B.getLength(), 14 )

class AlignatorDPWrapTestCase( AlignatorTestCase ):

    def setUp( self ):
        AlignatorTestCase.setUp( self )
        # penalize for gaps at the ends. This will force the residue to be aligned
        self.mAlignator = makeAlignatorDPFull( ALIGNMENT_WRAP, -10.0, -1.0 )

    def checkAlignment( self, row, col ):

        AlignatorTestCase.checkAlignment( self, row, col )

        if row == col or (row,col) in self.mSames:
            self.assertEqual( self.mAlignataA2B.getLength(), self.mSeqs[row].getLength() )            
        else:
            self.assertEqual( self.mAlignataA2B.getLength(), 14 )
            
def suite():
    suite = unittest.TestSuite()
    suite.addTest(AlignatorDPGlobalWithEndGapsPenaltiesTestCase)
    suite.addTest(AlignatorDPGlobalNoEndGapsPenaltiesTestCase)
    suite.addTest(AlignatorDPLocalTestCase)        
    return suite

if __name__ == "__main__":
    unittest.main()


