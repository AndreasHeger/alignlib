/*
  alignlib - a library for aligning protein sequences

  $Id: test_HelpersAlignment.cpp,v 1.4 2004/06/02 12:11:38 aheger Exp $

  Copyright (C) 2004 Andreas Heger
  
  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License
  as published by the Free Software Foundation; either version 2
  of the License, or (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
*/

/** Test alignata objects
*/

#if HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>

#include "alignlib.h"

using namespace std;
using namespace alignlib;

#define BOOST_TEST_MODULE
#include <boost/test/included/unit_test.hpp>
using boost::unit_test::test_suite;

BOOST_AUTO_TEST_CASE( test1 )
{
  const HEncoder translator = getDefaultEncoder();
  {
	  const HSubstitutionMatrix matrix = getDefaultSubstitutionMatrix();
    
	  std::string alphabet = translator->getAlphabet();
  
	  HAlignandum seq1 = makeSequence( alphabet );
    
	  for (int x = 0; x < seq1->getLength(); ++x)
	  {
		  Residue r = seq1->asResidue(x);
		  if (r == translator->getMaskCode() )
			  continue;
		  assert( matrix->getValue(r,r) > 0);
	  }
  }
}

// automatic mapping to 20 residue letter alphabet
BOOST_AUTO_TEST_CASE( test2 )
{
	const HEncoder t1(getEncoder( Protein20 ));
	const HEncoder t2(getEncoder( Protein23 ));	
	const HSubstitutionMatrix matrix1( makeSubstitutionMatrixBlosum62(t1) );
	const HSubstitutionMatrix matrix2( makeSubstitutionMatrixBlosum62() );
	BOOST_CHECK_EQUAL( matrix1->getNumRows(), t1->getAlphabetSize() );	
	BOOST_CHECK_EQUAL( matrix1->getNumCols(), t1->getAlphabetSize() );
	for (Residue x = 0; x < matrix1->getNumRows(); ++x )
		for (Residue y = 0; y < matrix1->getNumCols(); ++y )
			BOOST_CHECK_EQUAL( 
					matrix1->getValue(x,y), 
					matrix2->getValue( 
							t2->encode(t1->decode(x)),
							t2->encode(t1->decode(y)) ) );			
}	


// automatic mapping to 23 residue letter alphabet (no effect mapping).
BOOST_AUTO_TEST_CASE( test3 )
{
        const HEncoder t(getEncoder( Protein23 ));
        const HSubstitutionMatrix matrix1( makeSubstitutionMatrixBlosum62( t) );
        const HSubstitutionMatrix matrix2( makeSubstitutionMatrixBlosum62() );
        BOOST_CHECK_EQUAL( matrix1->getNumRows(), matrix2->getNumRows() );
        BOOST_CHECK_EQUAL( matrix1->getNumCols(), matrix2->getNumCols() );
        for (int x = 0; x < matrix1->getNumRows(); ++x )
                for (int y = 0; y < matrix1->getNumCols(); ++y )
                        BOOST_CHECK_EQUAL( matrix1->getValue(x,y), matrix2->getValue(x,y) );
}

// load substitution matrix from a file and compare to built-in
// This test only uses the twenty residue letter alphabet, because
// the blosum matrices are only defined for this alphabet and published
// matrices differ between scores to B/Z
BOOST_AUTO_TEST_CASE( test4 )
{
	const HEncoder t(getEncoder( Protein20 ));	  
	std::ifstream infile( "data/BLOSUM62");
	const HSubstitutionMatrix matrix_load( loadSubstitutionMatrix( infile, t ));
	const HSubstitutionMatrix matrix_ref( makeSubstitutionMatrixBlosum62(t) );
	BOOST_CHECK_EQUAL( matrix_load->getNumRows(), matrix_ref->getNumRows() );	
	BOOST_CHECK_EQUAL( matrix_load->getNumCols(), matrix_ref->getNumCols() );		
	for (int x = 0; x < matrix_load->getNumRows(); ++x )
		for (int y = 0; y < matrix_load->getNumCols(); ++y )
		{
			if ( matrix_load->getValue(x,y) != matrix_ref->getValue(x,y) )
				std::cout << "mismatch of scores: " << int(x) << " " << int(y) << " " << t->decode(x) << " " << t->decode(y) << std::endl;
			BOOST_CHECK_EQUAL( matrix_load->getValue(x,y), matrix_ref->getValue(x,y) );
		}
}	

// load substitution matrix from a file
BOOST_AUTO_TEST_CASE( test5 )
{
	const HEncoder t(getEncoder( Protein23 ));	  
	std::ifstream infile( "data/PAM30");
	const HSubstitutionMatrix matrix_load( loadSubstitutionMatrix( infile, t ));
	const HSubstitutionMatrix matrix_ref( makeSubstitutionMatrixPam30(t) );
	BOOST_CHECK_EQUAL( matrix_load->getNumRows(), matrix_ref->getNumRows() );	
	BOOST_CHECK_EQUAL( matrix_load->getNumCols(), matrix_ref->getNumCols() );		
	for (int x = 0; x < matrix_load->getNumRows(); ++x )
		for (int y = 0; y < matrix_load->getNumCols(); ++y )
		{
			if ( matrix_load->getValue(x,y) != matrix_ref->getValue(x,y) )
				std::cout << "mismatch of scores: " << int(x) << " " << int(y) << " " << t->decode(x) << " " << t->decode(y) << std::endl;
			BOOST_CHECK_EQUAL( matrix_load->getValue(x,y), matrix_ref->getValue(x,y) );
		}	
}	

// automatic mapping to 23 residue letter alphabet (no effect mapping).
BOOST_AUTO_TEST_CASE( test6 )
{
        const HEncoder t(getEncoder( Protein23 ));
        const HSubstitutionMatrix matrix1( makeSubstitutionMatrixBackTranslation( 1, -1, 0.5, t ) );
        std::cout << *matrix1 << std::endl;
        for (int x = 0; x < matrix1->getNumRows(); ++x )
                for (int y = 0; y < matrix1->getNumCols(); ++y )
                        BOOST_CHECK_EQUAL( matrix1->getValue(x,y), matrix1->getValue(y,x) );
        
}





