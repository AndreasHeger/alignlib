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
#include <vector>

#include "alignlib.h"
#include "MultAlignment.h"
#include "HelpersMultAlignment.h"

#define BOOST_TEST_MODULE
#include <boost/test/included/unit_test.hpp>
using boost::unit_test::test_suite;

using namespace std;
using namespace alignlib;

struct adder : public std::unary_function<double, void>
{
	adder() : sum(0) {}
	double sum;
	void operator()(double x) { sum += x; }
};

void testProfile( const HProfile & profile )
{
	// result should be
	// 0: no counts
	// 1: 3 * A
    // 2: 3 * A
	// 3: 3 * A
	// 4: 2 * A
	// 5: A + C + D
	// 6: 2 * C + D
	// 7: 2 * C + D
	// 8: C + D

	HCountMatrix matrix(profile->getCountMatrix());
	int sums[9]   = {0,3,3,3,2,3,3,3,2};
	int sums_A[9] = {0,3,3,3,2,1,0,0,0};
	int sums_C[9] = {0,0,0,0,0,1,2,2,1};
	int sums_D[9] = {0,0,0,0,0,1,1,1,1};
	int s = 9;

	BOOST_CHECK_EQUAL( matrix->getNumRows(), s);
	int size=matrix->getNumCols();
	HEncoder e = getDefaultEncoder();
	for( int x = 0; x < 9; ++x)
	{
		BOOST_CHECK_EQUAL( for_each( (*matrix)[x], (*matrix)[x]+size, adder()).sum, sums[x] );
		BOOST_CHECK_EQUAL( (*matrix)[x][e->encode('A')], sums_A[x] );
		BOOST_CHECK_EQUAL( (*matrix)[x][e->encode('C')], sums_C[x] );
		BOOST_CHECK_EQUAL( (*matrix)[x][e->encode('D')], sums_D[x] );
	}
}

// test creation of profile from two alignandum objects
BOOST_AUTO_TEST_CASE( test_makeProfile1a )
{
	HAlignandum a(makeSequence("AAAAACCCC"));
	HAlignandum b(makeProfile("AAAACCCCAAAADDDD", 2));
	HAlignment map_a2mali(makeAlignmentVector());
	map_a2mali->addDiagonal( 0,3,1);
	map_a2mali->addDiagonal( 4,7,1);
	HAlignment map_b2mali(makeAlignmentVector());
	map_b2mali->addDiagonal( 0,8,1);
	HProfile profile(toProfile( makeProfile( a, map_a2mali, b, map_b2mali)));
	testProfile( profile );
}

// test creation of profile from two alignandum objects
BOOST_AUTO_TEST_CASE( test_makeProfile1b )
{
	HAlignandum a(makeSequence("AAAAACCCC"));
	HAlignandum b(makeProfile("AAAACCCCAAAADDDD", 2));
	HAlignment map_a2mali(makeAlignmentVector());
	map_a2mali->addDiagonal( 0,3,1);
	map_a2mali->addDiagonal( 4,7,1);
	HAlignment map_b2mali(makeAlignmentVector());
	map_b2mali->addDiagonal( 0,8,1);
	HProfile profile(toProfile( makeProfile( b, map_b2mali, a, map_a2mali)));
	testProfile( profile );
}

// test creation of profile from a MultAlignment
BOOST_AUTO_TEST_CASE( test_makeProfile2a )
{
	HAlignandum a(makeSequence("AAAAACCCC"));
	HAlignandum b(makeProfile("AAAACCCCAAAADDDD", 2));
	HAlignment map_a2mali(makeAlignmentVector());
	map_a2mali->addDiagonal( 0,3,1);
	map_a2mali->addDiagonal( 4,7,1);
	map_a2mali->switchRowCol();
	HAlignment map_b2mali(makeAlignmentVector());
	map_b2mali->addDiagonal( 0,8,1);
	map_b2mali->switchRowCol();
	HMultAlignment mali(makeMultAlignment());
	mali->add( map_a2mali );
	mali->add( map_b2mali );
	HAlignandumVector seqs(new AlignandumVector());
	seqs->push_back( a );
	seqs->push_back( b );

	HProfile profile(toProfile( makeProfile( mali, seqs)));
	testProfile( profile );
}

// test creation of profile from a MultAlignment
BOOST_AUTO_TEST_CASE( test_makeProfile2b )
{
	HAlignandum a(makeSequence("AAAAACCCC"));
	HAlignandum b(makeProfile("AAAACCCCAAAADDDD", 2));
	HAlignment map_a2mali(makeAlignmentVector());
	map_a2mali->addDiagonal( 0,3,1);
	map_a2mali->addDiagonal( 4,7,1);
	map_a2mali->switchRowCol();
	HAlignment map_b2mali(makeAlignmentVector());
	map_b2mali->addDiagonal( 0,8,1);
	map_b2mali->switchRowCol();
	HMultAlignment mali(makeMultAlignment());
	mali->add( map_b2mali );
	mali->add( map_a2mali );
	HAlignandumVector seqs(new AlignandumVector());
	seqs->push_back( b );
	seqs->push_back( a );

	HProfile profile(toProfile( makeProfile( mali, seqs)));
	testProfile( profile );
}










