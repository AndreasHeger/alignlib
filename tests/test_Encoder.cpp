/*
  alignlib - a library for aligning protein sequences

  $Id: test_MultipleAlignment.cpp,v 1.6 2004/06/02 12:11:38 aheger Exp $

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

/** Test the Encoder object
 */

#if HAVE_CONFIG_H
#include <config.h>
#endif

#include <cstdio>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cassert>

#define BOOST_TEST_MODULE
#include <boost/test/included/unit_test.hpp>
using boost::unit_test::test_suite;

#include "Encoder.h"
#include "HelpersEncoder.h"

using namespace std;
using namespace alignlib;

void testEncoder( const HEncoder & translator,
		const std::string & alphabet,
		const std::string & gap_chars,
		const std::string & mask_chars,
		int alphabet_size )
{
	// check alphabet size (number of chars + 1 for mask value)
	BOOST_CHECK_EQUAL( alphabet_size, translator->getAlphabetSize() );

	for (int x = 0; x < alphabet.size(); ++x)
	{
		BOOST_CHECK_EQUAL( x, translator->encode( alphabet[x] ) );
		// check if upper/lower-case is ignored
		BOOST_CHECK_EQUAL( translator->encode( toupper(alphabet[x]) ), translator->encode( tolower( alphabet[x] ) ) );
	}

	// check mask codes
	for (int x = 0; x < mask_chars.size(); ++ x)
		BOOST_CHECK_EQUAL( translator->encode(mask_chars[x]), translator->getMaskCode() );

	// check gap codes
	for (int x = 0; x < gap_chars.size(); ++ x)
		BOOST_CHECK_EQUAL( translator->encode(gap_chars[x]), translator->getGapCode() );

	// check saving/loading from stream.
	{
		ofstream file("test_Encoder.tmp", ios::binary);
		translator->save( file );
		translator->save( file );
		file.close();
	}

	{
		ifstream file("test_Encoder.tmp", ios::binary) ;

		HEncoder b;
		int n = 0;
		while ( file.peek() != EOF )
		{
			b = loadEncoder( file ) ;
			if (b->getAlphabetType() != User )
				BOOST_CHECK_EQUAL( translator, b );
			++n;
		}
		BOOST_CHECK_EQUAL( n, 2 );
		std::remove( "test_Encoder.tmp" );
	}

	{
		HEncoder n( getEncoder( Protein20 ));
		HResidueVector map_new2old( n->map( translator ) );
		for (Residue x = 0; x < map_new2old->size(); ++x)
			{}
			// std::cout << (int)x << "\t" << (int)(*map_new2old)[x] << std::endl;
	}

	{
		HEncoder n( getEncoder( Protein20 ));
		HResidueVector map_new2old( translator->map( n ) );
		for (Residue x = 0; x < map_new2old->size(); ++x)
		{}
		// std::cout << (int)x << "\t" << (int)(*map_new2old)[x] << std::endl;
	}

}

BOOST_AUTO_TEST_CASE( testProtein23 )
{
	testEncoder( getEncoder( Protein23 ), "ABCDEFGHIKLMNPQRSTVWXYZ", "-.", "X", 23 );
}

BOOST_AUTO_TEST_CASE( testProtein20 )
{
	testEncoder( getEncoder( Protein20 ), "ACDEFGHIKLMNPQRSTVWY", "-.", "X", 21 );
}

BOOST_AUTO_TEST_CASE( testDNA4 )
{
	testEncoder( getEncoder( DNA4 ), "ACGT", "-.", "N", 5 );
}

