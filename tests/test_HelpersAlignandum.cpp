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

#include <time.h>

#include "alignlib.h"
#include "alignlib_fwd.h"

using namespace std;
using namespace alignlib;

#define BOOST_TEST_MODULE
#include <boost/test/included/unit_test.hpp>
using boost::unit_test::test_suite;

BOOST_AUTO_TEST_CASE( test_makeSequenceFromFasta )
{
	std::vector<std::string>in;
	std::vector<std::string>out;

	in.push_back( "ACT");
	out.push_back( "ACT");

	in.push_back( "ACT\nACT");
	out.push_back( "ACTACT");

	in.push_back( "A C T\n A C T");
	out.push_back( "ACTACT");

	{
		ofstream outfile("test_HelpersAlignandum.tmp");
		for (int x = 0; x < in.size(); ++x)
			outfile << ">" << x << "\n" << in[x] << std::endl;
		outfile.close();
	}

	{
		ifstream file("test_HelpersAlignandum.tmp") ;
		std::string description;
		HAlignandum b;
		int x = 0;
		while ( file.peek() != EOF )
		{
			b = makeSequenceFromFasta( file, description );
			BOOST_CHECK_EQUAL( atoi(description.c_str()), x);
			BOOST_CHECK_EQUAL( out[x], b->asString() );
			++x;
		}
		assert( x == in.size() );
	}
}

