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

#define BOOST_TEST_MODULE
#include <boost/test/included/unit_test.hpp>
using boost::unit_test::test_suite;

using namespace std;
using namespace alignlib;

HMultAlignment buildAlignment()
{
	HMultAlignment a(makeMultAlignment());
	HAlignment ali(makeAlignmentVector());
	a->add( ali );
	for (int x = 0; x < 10; ++x)
	{
		ali->clear();
		ali->addDiagonal( x, x+10, -x);
		a->add( ali );
	}
	return a;
}


void testWriteRead( const std::auto_ptr<MultAlignmentFormat> & format,
					const HMultAlignment & old_mali )
{
	stringstream ss;

	// write
	ss << *format;
	ss << std::ends;

	// read
	MultAlignmentFormatPlain new_format;
	ss >> new_format;
	HMultAlignment new_mali(old_mali->getNew());
	new_format.copy( new_mali, makeAlignmentVector() );

	// compare
	BOOST_CHECK_EQUAL( format->mData.size(), new_format.mData.size() );

	for (int x = 0; x < format->mData.size(); ++x)
		BOOST_CHECK_EQUAL(
				format->mData[x]->getString(),
				new_format.mData[x]->getString());

	BOOST_CHECK_EQUAL( checkMultAlignmentIdentity( old_mali, new_mali ), true );
}

BOOST_AUTO_TEST_CASE( test_MultAlignmentFormatPlain )
{
	HMultAlignment mali(buildAlignment());
	HStringVector sequences(new StringVector());
	for (int x = 0; x < mali->getNumSequences(); ++x)
	{
		sequences->push_back( "XXXXXXXXXXX" );
	}

	std::auto_ptr<MultAlignmentFormat>f(new MultAlignmentFormatPlain( mali, sequences));
	testWriteRead( f, mali );
}



