# alignlib - a library for aligning protein sequences
# 
# $Id$
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

class RegularizorCase( unittest.TestCase ):

    def setUp( self ):
        
        self.mFactory = makeRegularizor

        self.mMali = makeMultipleAlignment()
        
        self.mMali.add( makeAlignatum("ACDEFGHIKLMNPQRSTVWY") )
        self.mMali.add( makeAlignatum("ACDEFGHIKLMNPQRSTVWY") )
        self.mMali.add( makeAlignatum("ACDEFGHIKLMNPQRSTVWY") )
        
    def testRegularizor( self ):
                            
        # getDefaultToolkit().setEncoder( getEncoder( Protein20) )
        getDefaultToolkit().setRegularizor( self.mFactory() )
        a = makeProfile(self.mMali)
        a.prepare()
        print a

class RegularizorDirichletCase( RegularizorCase ):

    def setUp( self ):
        RegularizorCase.setUp( self )
        self.mFactory = makeRegularizorDirichletPrecomputed

"""
class WeightorHenikoffCase( WeightorCase ):

    def setUp( self ):
        self.mWeightorFactory = makeWeightorHenikoff

class WeightorHenikoffKimmenCase( WeightorCase ):

    def setUp( self ):
        self.mWeightorFactory = makeWeightorHenikoffKimmen
"""     

class TatusovCase( unittest.TestCase ):

    def testPsiblast( self ):

        matrix = makeSubstitutionMatrixBlosum62()

        bg = (0.078047, 0.053640, 0.062949, 0.038556, 0.038556,                                                                                                                                                              
                0.073772, 0.021992, 0.051420, 0.057438, 0.090191,                                                                                                                                  
                0.022425, 0.044873, 0.052028, 0.042644, 0.051295,                                                                                                                                  
                0.071198, 0.058413, 0.064409, 0.013298, 0.032165 )

        b = FrequencyVector()
        b.extend( bg )
        regularizor = makeRegularizorTatusov(
            matrix,
            b,
            "ACDEFGHIKLMNPQRSTVWY",                                                                                                                                                                                              
            10, 
            0.3176 )

    def testDNA( self ):

        matrix = makeSubstitutionMatrixDNA4()

        bg = (0.4, 0.3, 0.2, 0.2)

        weight_matrix = MatrixDouble( 10, 5, 0.0 )

        for x in range(0, 10):
            weight_matrix.setValue(x, x % 5, 1.0 )

        printMatrix( sys.stdout, weight_matrix )

        frequency_matrix = MatrixDouble( 10, 5, 0.0 )

        b = FrequencyVector()
        b.extend( bg )
        regularizor = makeRegularizorTatusov(
            matrix,
            b,
            "ACGT",                                                                                                                                                                                              
            10.0, 
            0.0 )

        regularizor.fillFrequencies( frequency_matrix,
                                     weight_matrix,
                                     getEncoder( DNA4) )
        
        printMatrix( sys.stdout, frequency_matrix )


def suite():
    suite = unittest.TestSuite()
    suite.addTest(RegularizorDirichletCase)
    suite.addTest(TatusovCase)
    return suite

if __name__ == "__main__":
    unittest.main()

