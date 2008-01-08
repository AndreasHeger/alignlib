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
#include "Translator.h"
#include "HelpersTranslator.h"
#include "HelpersAlignandum.h"
#include "HelpersSubstitutionMatrix.h"
#include "HelpersProfile.h"

#define BOOST_TEST_MODULE
#include <boost/test/included/unit_test.hpp>
using boost::unit_test::test_suite;

using namespace std;
using namespace alignlib;

std::string ref_protein20 = "ACDEFGHIKLMNPQRSTVWY"; 
std::string ref_protein20x3 = ref_protein20 + ref_protein20 + ref_protein20; 

// most profilers only work with the twenty amino acid alphabet
// setDefaultTranslator( getTranslator(Protein20) );

void test_GenericRegularizor( const HRegularizor & r )
{	
	std::auto_ptr<Alignandum>a(makeProfile( ref_protein20x3, 3,
			getDefaultTranslator(), 
			getDefaultWeightor(), 
			makeRegularizor(), 
			getDefaultLogOddor() ));	
	a->prepare();
	std::cout << *a << std::endl;
}

BOOST_AUTO_TEST_CASE( test_Regularizor )
{
	HRegularizor r( makeRegularizor());
	test_GenericRegularizor( r );
}

BOOST_AUTO_TEST_CASE( test_RegularizorPsiblast )
{
	HRegularizor r( makeRegularizorPsiblast());
	test_GenericRegularizor( r );
}








