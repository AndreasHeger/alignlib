/*
  alignlib - a library for aligning protein sequences

  $Id: test_HelpersAlignment.cpp,v 1.5 2004/10/14 23:34:09 aheger Exp $

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
#include <sstream>
#include <cstdio>
#include <time.h> 

#include "alignlib.h"
#include "alignlib_fwd.h"
#include "AlignlibIndex.h"

using namespace std;
using namespace alignlib;

#define BOOST_TEST_MODULE
#include <boost/test/included/unit_test.hpp>
using boost::unit_test::test_suite;

typedef Index< int, RecorderTable< int, 0 > > IntIndex; 
typedef Index< int, RecorderTable< int, 1 > > IntIndex2; 
typedef Index< std::string, RecorderTable< std::string, 0 > > StringIndex;

BOOST_AUTO_TEST_CASE( test_IndexCreateColumn1 )
{
	FILE * file = openFileForRead( "./data/test_index.data");
	IntIndex index( file );
	BOOST_CHECK_EQUAL( 3, index.size() );
	std::fclose(file);
}

BOOST_AUTO_TEST_CASE( test_IndexCreateColumn2 )
{
	FILE * file = openFileForRead( "./data/test_index.data");
	IntIndex2 index( file );
	BOOST_CHECK_EQUAL( 6, index.size() );
	std::fclose(file);
}

BOOST_AUTO_TEST_CASE( test_IndexPosition )
{
	FILE * file = openFileForRead( "./data/test_index.data");
	IntIndex index( file );
	char buffer[1000];
	int r = 0;
	for( int x = 1; x < 4; ++x)
	{
		index.goTo( x );
		std::fgets( buffer, 1000, file );
		std::istringstream i( buffer );
		i >> r;
		BOOST_CHECK_EQUAL( r, x );		
	}
	std::fclose(file);
}
