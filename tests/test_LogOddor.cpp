/*
  alignlib - a library for aligning protein sequences

  $Id: test_Alignandum.cpp,v 1.3 2004/06/02 12:11:38 aheger Exp $

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
#include <fstream>
#include <cassert>
#include <time.h> 

#include "alignlib.h"
#include "alignlib_fwd.h"

#include "Alignandum.h"
#include "Encoder.h"
#include "Matrix.h"
#include "HelpersEncoder.h"
#include "HelpersAlignandum.h"
#include "HelpersSubstitutionMatrix.h"

#define BOOST_TEST_MODULE
#include <boost/test/included/unit_test.hpp>
using boost::unit_test::test_suite;

using namespace std;
using namespace alignlib;

std::string ref_protein20 = "ACDEFGHIKLMNPQRSTVWY";
std::string ref_protein20A = "AAAAAAAAAAAAAAAAAAAA";

std::string ref_protein20x3 = ref_protein20 + ref_protein20 + ref_protein20 + ref_protein20A; 

void test_GenericLogOddor( HLogOddor & logoddor )
{
	HAlignandum a(makeProfile( ref_protein20x3, 4,
			getDefaultEncoder(), 
			getDefaultWeightor(), 
			getDefaultRegularizor(), 
			logoddor ));
	
	a->prepare();	
	std::cout << *a << std::endl;
}

BOOST_AUTO_TEST_CASE( test_LogOddor )
{
	HLogOddor l = makeLogOddor();	
	test_GenericLogOddor( l );
}

BOOST_AUTO_TEST_CASE( test_LogOddorUniform )
{
	HLogOddor l = makeLogOddorUniform();
	test_GenericLogOddor( l );
}

BOOST_AUTO_TEST_CASE( test_LogOddorGribskov )
{
	setDefaultEncoder( getEncoder( Protein23 ) );
	HSubstitutionMatrix m( makeSubstitutionMatrixBlosum62() );
	HLogOddor l = makeLogOddorGribskov( m );
	test_GenericLogOddor( l ); 
}

BOOST_AUTO_TEST_CASE( test_LogOddorDirichlet20 )
{
	setDefaultEncoder( getEncoder( Protein20 ) );
	HLogOddor l = makeLogOddorDirichlet();
	test_GenericLogOddor( l );
}

BOOST_AUTO_TEST_CASE( test_LogOddorDirichlet23 )
{
	setDefaultEncoder( getEncoder( Protein23 ) );
	HLogOddor l = makeLogOddorDirichlet();
	test_GenericLogOddor( l );
}






