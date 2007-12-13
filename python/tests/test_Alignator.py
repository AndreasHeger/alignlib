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
import sys

class AlignatorTestCase( unittest.TestCase ):

    def __init__(self, *args, **kwargs):
        self.mReferenceSequence1 = None
        self.mReferenceSequence2 = None
        unittest.TestCase.__init__(self, *args, **kwargs )

    def setUp( self ):

        if not self.mReferenceSequence1:
            self.mReferenceSequence1 = "AAAAAAACCCCAAAAAAA"
        if not self.mReferenceSequence2:
            self.mReferenceSequence2 = "WWAAAAAWWWWAAAAAWW"
                     
        self.mAlignmentA2B = makeAlignmentVector()
        self.mAlignmentB2A = makeAlignmentVector()        
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
            self.mAlignmentA2B.clear()
            self.mAlignmentB2A.clear()            
            self.mAlignator.align( self.mSeqs[row], self.mSeqs[col], self.mAlignmentA2B )
            self.mAlignator.align( self.mSeqs[col], self.mSeqs[row], self.mAlignmentB2A )            
            self.checkAlignment( row, col )

    def checkAlignment( self, row, col ):
        """general sanity checks."""

        self.assertEqual( self.mAlignmentA2B.getLength(), self.mAlignmentB2A.getLength() )
        self.assertEqual( self.mAlignmentA2B.getScore(), self.mAlignmentB2A.getScore() )
        self.assertEqual( self.mAlignmentA2B.getNumGaps(), self.mAlignmentB2A.getNumGaps() )
        self.assertEqual( self.mAlignmentA2B.getRowFrom(), self.mAlignmentB2A.getColFrom() )
        self.assertEqual( self.mAlignmentA2B.getRowTo(), self.mAlignmentB2A.getColTo() )        
        
class AlignatorDPGlobalWithEndGapsPenaltiesTestCase( AlignatorTestCase ):

    def setUp( self ):
        AlignatorTestCase.setUp( self )
        # penalize for gaps at the ends. This will force the residue to be aligned
        self.mAlignator = makeAlignatorDPFull( ALIGNMENT_GLOBAL, -10.0, -1.0, True, True )

    def checkAlignment( self, row, col ):

        AlignatorTestCase.checkAlignment( self, row, col )
        self.assertEqual( self.mAlignmentA2B.getLength(), self.mSeqs[row].getLength() )            

class AlignatorDPGlobalNoEndGapsPenaltiesTestCase( AlignatorTestCase ):

    def setUp( self ):
        AlignatorTestCase.setUp( self )
        # do not penalize for gaps at the ends (default). This will allow mismatching residues
        # at the ends to be unaligned
        self.mAlignator = makeAlignatorDPFull( ALIGNMENT_GLOBAL, -10.0, -1.0, False, False )

    def checkAlignment( self, row, col ):

        AlignatorTestCase.checkAlignment( self, row, col )
        
        if row == col or (row,col) in self.mSames:
            self.assertEqual( self.mAlignmentA2B.getLength(), self.mSeqs[row].getLength() )            
        else:
            self.assertEqual( self.mAlignmentA2B.getLength(), self.mSeqs[row].getLength() - 4 )

class AlignatorDPLocalTestCase( AlignatorTestCase ):

    def setUp( self ):
        AlignatorTestCase.setUp( self )
        # penalize for gaps at the ends. This will force the residue to be aligned
        self.mAlignator = makeAlignatorDPFull( ALIGNMENT_LOCAL, -10.0, -1.0 )

    def checkAlignment( self, row, col ):

        AlignatorTestCase.checkAlignment( self, row, col )

        if row == col or (row,col) in self.mSames:
            self.assertEqual( self.mAlignmentA2B.getLength(), self.mSeqs[row].getLength() )            
        else:
            self.assertEqual( self.mAlignmentA2B.getLength(), 14 )

class AlignatorDPWrapTestCase( AlignatorTestCase ):

    def setUp( self ):
        AlignatorTestCase.setUp( self )
        # penalize for gaps at the ends. This will force the residue to be aligned
        self.mAlignator = makeAlignatorDPFull( ALIGNMENT_WRAP, -10.0, -1.0 )

    def checkAlignment( self, row, col ):

        AlignatorTestCase.checkAlignment( self, row, col )

        if row == col or (row,col) in self.mSames:
            self.assertEqual( self.mAlignmentA2B.getLength(), self.mSeqs[row].getLength() )            
        else:
            self.assertEqual( self.mAlignmentA2B.getLength(), 14 )
            

class AlignatorIterativeTestCase( AlignatorTestCase ):

    def setUp( self ):

        self.mReferenceSequence1 = "AAACCCCCCCCGGGCCCCCCCQQQCCCCCCCWWWCCCCCCCFFF"
        self.mReferenceSequence2 = "AAAKKKKKKKkGGGKKKKKKKQQQKKKKKKKWWWKKKKKKKFFF"        

        self.mReferenceSequence1 = "AAACCCCCCCCGGGCCCCCCCQQQCCCCCCCAAACCCCCCCFFF"
        self.mReferenceSequence2 = "AAAKKKKKKKkGGGKKKKKKKQQQKKKKKKKAAAKKKKKKKFFF"        

        AlignatorTestCase.setUp( self )
        
        alignator = makeAlignatorDPFull( ALIGNMENT_LOCAL, -10.0, -2.0 )
        
        self.mAlignator = makeAlignatorIterative( alignator, 1.0 )
        
    def checkAlignment( self, row, col ):

        if self.mAlignmentA2B.getScore() != self.mAlignmentB2A.getScore():
            print str(self.mAlignmentA2B)
            print str(self.mAlignmentB2A)
            print row, col, self.mAlignmentA2B.getScore(), self.mAlignmentB2A.getScore()

        AlignatorTestCase.checkAlignment( self, row, col )
        
        if row == col or (row,col) in self.mSames:
            self.assertEqual( self.mAlignmentA2B.getLength(), self.mSeqs[row].getLength() )            
            self.assertEqual( self.mAlignmentA2B.getLength(), self.mSeqs[col].getLength() )
            self.assertEqual( self.mAlignmentA2B.getNumGaps(), 0 )
        
def suite():
    suite = unittest.TestSuite()
    suite.addTest(AlignatorDPGlobalWithEndGapsPenaltiesTestCase)
    suite.addTest(AlignatorDPGlobalNoEndGapsPenaltiesTestCase)
    suite.addTest(AlignatorDPLocalTestCase)
    suite.addTest(AlignatorIterativeTestCase)   
    return suite

if __name__ == "__main__":
    unittest.main()


