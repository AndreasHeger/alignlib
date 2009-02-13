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
import unittest, sys, os
from alignlib import *

class MultAlignmentTestCase( unittest.TestCase ):

    mReferenceSequence = "0123456789"
    mNumSequences = 3
        
    def setUp( self ):
        self.mAlignandum = makeSequence( self.mReferenceSequence )
        self.mContainer = makeAlignmentVector()
        
    def constructMali(self):
        mali = makeMultAlignment()
        ali = self.mContainer.getNew()
        ali.addDiagonal( 0,3,+2 );
        ali.addDiagonal( 3,6,+4 );
        mali.add( ali );

        ali = self.mContainer.getNew()
        ali.addDiagonal( 0,1,+1 );
        ali.addDiagonal( 1,6,+3 );
        mali.add( ali );
        mali.add( ali );
        
        seqs = StringVector()
        for x in range( self.mNumSequences):
            seqs.append( self.mReferenceSequence )

        return mali, seqs
    
    def testBuild(self):
        mali, seqs = self.constructMali()
        self.assertEqual( mali.getNumSequences(), len(seqs) )
        self.assertEqual( mali.getLength(), 6 )

    def testExpandSimple(self):
        """expand mali without sequences."""
        mali, seqs = self.constructMali()
        mali.expand( AlignandumVector() )
        format = MultAlignmentFormatPlain( mali, seqs )
        result = [ x.split("\t") for x in str(format).split("\n") ]
        self.assertEqual( result[0], ["2", "2----3456789", "10" ] )
        self.assertEqual( result[1], ["1", "123--45--678", "10" ] )
        self.assertEqual( result[2], ["1", "1--2345--678", "10" ] )

    def testExpandFull(self):
        """expand mali with sequences."""
        mali, seqs = self.constructMali()
        v = AlignandumVector()
        for x in seqs: v.append( makeSequence(x) )
        mali.expand( v )
        format = MultAlignmentFormatPlain( mali, seqs )
        result = [ x.split("\t") for x in str(format).split("\n") ]
        self.assertEqual( result[0], ["0", "01--2----3456789--", "10" ] )
        self.assertEqual( result[1], ["0", "--0-123--45--6789-", "10" ] )
        self.assertEqual( result[2], ["0", "---01--2345--678-9", "10" ] )

class MultAlignmentBlocksTestCase( MultAlignmentTestCase ):
    def setUp( self ):
        MultAlignmentTestCase.setUp( self )
        self.mContainer = makeAlignmentBlocks()

class MultAlignmentSetTestCase( MultAlignmentTestCase ):
    def setUp( self ):
        MultAlignmentTestCase.setUp( self )
        self.mContainer = makeAlignmentSet()

class MultAlignmentHashTestCase( MultAlignmentTestCase ):
    def setUp( self ):
        MultAlignmentTestCase.setUp( self )
        self.mContainer = makeAlignmentHash()

class MultAlignmentSetColTestCase( MultAlignmentTestCase ):
    def setUp( self ):
        MultAlignmentTestCase.setUp( self )
        self.mContainer = makeAlignmentSetCol()

class MultAlignmentHashDiagonalTestCase( MultAlignmentTestCase ):
    def setUp( self ):
        MultAlignmentTestCase.setUp( self )
        self.mContainer = makeAlignmentHashDiagonal()

def suite():
    suite = unittest.TestSuite()
    suite.addTest(MultAlignmentTestCase)
    return suite

if __name__ == "__main__":
    unittest.main()


        





