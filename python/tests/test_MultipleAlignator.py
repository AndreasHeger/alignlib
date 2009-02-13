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

class MultipleAlignatorTestCase( unittest.TestCase ):

    def __init__(self, *args, **kwargs):
        unittest.TestCase.__init__(self, *args, **kwargs )

        setDefaultSubstitutionMatrix( makeSubstitutionMatrix( getDefaultEncoder().getAlphabetSize(), 1, -1000) )
        self.mSequences = StringVector()
        self.mSequences.append( "EEEEAAAADDDDMMMMEEEE");
        self.mSequences.append( "FFFFAAAACCCCMMMMFFFF");
        self.mSequences.append( "GGGGAAAAKKKKMMMMGGGG");
        self.mResult = makeMultAlignment()
        self.mAlignator = makeAlignatorDPFull( ALIGNMENT_LOCAL, 0, 0 )  
        self.mMultipleAlignator = None
        
    def testAlignment(self):
        if self.mMultipleAlignator:
            self.mMultipleAlignator.align( self.mResult, self.mSequences )   
            f = MultAlignmentFormatPlain( self.mResult, self.mSequences)
            for x in range( 0, len(self.mSequences)):                        
                self.assertEqual( f.mData[x].getString(), self.mResults[x] )
                
class MultipleAlignatorSimpleTestCase( MultipleAlignatorTestCase ):

    def __init__(self, *args, **kwargs):
        MultipleAlignatorTestCase.__init__(self, *args, **kwargs )
        self.mMultipleAlignator = makeMultipleAlignatorSimple( self.mAlignator ) 
        self.mResults = ["--------EEEEAAAADDDD--------MMMMEEEE--------",
                         "----FFFF----AAAA----CCCC----MMMM----FFFF----",
                         "GGGG--------AAAA--------KKKKMMMM--------GGGG" ]
        
class MultipleAlignatorPileupTestCase( MultipleAlignatorTestCase ):

    def __init__(self, *args, **kwargs):
        MultipleAlignatorTestCase.__init__(self, *args, **kwargs )
                
        self.mMultipleAlignator = makeMultipleAlignatorPileup( self.mAlignator ) 
        self.mResults = ["EEEEAAAADDDDMMMMEEEE",
                         "----AAAA----MMMM----",
                         "----AAAA----MMMM----", ]
                        
def suite():
    suite = unittest.TestSuite()
    suite.addTest( MultipleAlignatorSimpleTestCase )
    suite.addTest( MultipleAlignatorPileupTestCase )
    return suite

if __name__ == "__main__":
    unittest.main()

