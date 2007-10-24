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

class AlignandumTestCase( unittest.TestCase ):

    mReferenceSequence = "AAACCCCAAA"

    def setUp( self ):
        
        self.mAlignandum = makeSequence( self.mReferenceSequence )
        
    def renderSeq( self ):
        return str(self.mAlignandum)

    def testAsString( self ):
        self.assertEqual( self.mReferenceSequence, self.mAlignandum.asString() )
        
    def testStr( self ):
        self.assertEqual( self.mReferenceSequence, self.renderSeq() )

    def testAsChar( self ):

        for x in range(len(self.mReferenceSequence)):
            self.assertEqual( self.mReferenceSequence[x], self.mAlignandum.asChar(x+1) )

    def testGetLength( self ):
        self.assertEqual( self.mAlignandum.getLength(), 10 )

    def testPrepare( self ):
        # a sequence is always prepared
        self.mAlignandum.prepare()
        self.assertEqual( self.mAlignandum.isPrepared(), True )
        self.mAlignandum.release()
        self.assertEqual( self.mAlignandum.isPrepared(), True )        

    def testUseSegment( self ):

        for start, end in ( (0,len(self.mReferenceSequence)),
                            (4,len(self.mReferenceSequence)) ):
            self.mAlignandum.useSegment( start+1, end)
            self.assertEqual( self.mReferenceSequence[start:end],
                              self.mAlignandum.asString() )

        self.mAlignandum.useFullLength()
        self.assertEqual( self.mReferenceSequence,
                          self.mAlignandum.asString() )
        
    def testMask( self ):
        self.mAlignandum.mask( 4, 7 )
        self.assertEqual( "AAAXXXXAAA", self.renderSeq() )

class ProfileTestCase( AlignandumTestCase ):

    def setUp( self ):
        n = len(self.mReferenceSequence)

        s = []
        for i in range (1, n + 1) :
            s.append( self.mReferenceSequence[0:i] + "-" * (n - i) )
        self.mAlignandum = makeProfile( "".join(s), n )

    def renderSeq( self ):
        return self.mAlignandum.asString()

    def testPrepare( self ):
        self.mAlignandum.prepare()
        self.assertEqual( self.mAlignandum.isPrepared(), True )
        self.mAlignandum.release()
        self.assertEqual( self.mAlignandum.isPrepared(), False )        

    def testMask( self ):
        ## TODO: make this test and the code conforming to use X as mask char
        self.mAlignandum.mask( 4, 7 )
        self.assertEqual( "AAA----AAA", self.renderSeq() )

def suite():
    suite = unittest.TestSuite()
    suite.addTest(AlignandumTestCase)
    suite.addTest(ProfileTestCase)    
    return suite

if __name__ == "__main__":
    unittest.main()


        





