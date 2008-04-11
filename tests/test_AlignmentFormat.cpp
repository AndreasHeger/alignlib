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

HAlignment buildAlignment()
{
	HAlignment a(makeAlignmentVector());
	addDiagonal2Alignment( a, 2, 4, 0);
	addDiagonal2Alignment( a, 6, 8, 0);
	return a;
}

void testWriteRead( const std::auto_ptr<AlignmentFormat> & format,
					const HAlignment & ref )
{
	// stringstream ss;
	// ss << *format;
		
	HAlignment n(ref->getNew());
	format->copy( n );
	BOOST_CHECK( checkAlignmentIdentity( ref, n ) );	
}

BOOST_AUTO_TEST_CASE( test_AlignmentExplicit )
{
	HAlignment ali(buildAlignment());
	HAlignandum s1(makeSequence( "AADDEEGGLL" ) );
	HAlignandum s2(makeSequence( "CCDDFFGGMM" ) );	
	
	std::auto_ptr<AlignmentFormat>f(new AlignmentFormatExplicit( ali, s1, s2 ));
	testWriteRead( f, ali );
}

BOOST_AUTO_TEST_CASE( test_AlignmentBlocks )
{
	HAlignment ali(buildAlignment());
	std::auto_ptr<AlignmentFormat>f(new AlignmentFormatBlocks( ali));	
	testWriteRead( f, ali );
}

BOOST_AUTO_TEST_CASE( test_AlignmentBlat )
{
	HAlignment ali(buildAlignment());
	std::auto_ptr<AlignmentFormat>f(new AlignmentFormatBlat( ali));	
	testWriteRead( f, ali );
}

BOOST_AUTO_TEST_CASE( test_AlignmentEmissions )
{
	HAlignment ali(buildAlignment());
	std::auto_ptr<AlignmentFormat>f(new AlignmentFormatEmissions( ali));	
	testWriteRead( f, ali );
}

BOOST_AUTO_TEST_CASE( test_AlignmentDiagonals )
{
	HAlignment ali(buildAlignment());	
	std::auto_ptr<AlignmentFormat>f(new AlignmentFormatDiagonals( ali));	
	testWriteRead( f, ali );
}


