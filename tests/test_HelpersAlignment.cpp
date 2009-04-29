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

BOOST_AUTO_TEST_CASE( test_fillAlignmentGaps1 )
{
	HAlignment ali( makeAlignmentVector() );
	addDiagonal2Alignment( ali, 0, 5, 0);
	addDiagonal2Alignment( ali, 6,10, 0);
	fillAlignmentGaps( ali, 3 );
	BOOST_CHECK_EQUAL( ali->getNumGaps(), 0);
	BOOST_CHECK_EQUAL( ali->getLength(), 10);
}

// test filling an empty alignment
BOOST_AUTO_TEST_CASE( test_fillAlignmentGaps2 )
{
	HAlignment ali( makeAlignmentVector() );
	fillAlignmentGaps( ali, 3 );
	BOOST_CHECK_EQUAL( ali->getNumGaps(), 0);
	BOOST_CHECK_EQUAL( ali->getLength(), 0);
}

BOOST_AUTO_TEST_CASE( test_hasAlignmentOverlap )
{
	HAlignment ali1( makeAlignmentVector() );
	ali1->addDiagonal( 0, 10, 20);

	for (int x = 0; x < 10; ++x)
	{
		HAlignment ali2( makeAlignmentVector() );
		ali2->addDiagonal( x, x+10, 20);
		BOOST_CHECK_EQUAL( hasAlignmentOverlap( ali1, ali2, RR, 1), true);
		BOOST_CHECK_EQUAL( hasAlignmentOverlap( ali1, ali2, CR, 1), false);
		BOOST_CHECK_EQUAL( hasAlignmentOverlap( ali1, ali2, RC, 1), false);
		BOOST_CHECK_EQUAL( hasAlignmentOverlap( ali1, ali2, CC, 1), true);
	}
}

BOOST_AUTO_TEST_CASE( test_getAlignmentOverlap )
{
	HAlignment ali1( makeAlignmentVector() );
	ali1->addDiagonal( 0, 10, 20);

	for (int x = 0; x < 10; ++x)
	{
		HAlignment ali2( makeAlignmentVector() );
		ali2->addDiagonal( x, x+10, 20);
		BOOST_CHECK_EQUAL( getAlignmentOverlap( ali1, ali2, RR ), 10 - x);
		BOOST_CHECK_EQUAL( getAlignmentOverlap( ali1, ali2, CR ), 0);
		BOOST_CHECK_EQUAL( getAlignmentOverlap( ali1, ali2, RC ), 0);
		BOOST_CHECK_EQUAL( getAlignmentOverlap( ali1, ali2, CC ), 10 - x);
	}
}

BOOST_AUTO_TEST_CASE( test_getAlignmentIdentity )
{
	HAlignment ali1( makeAlignmentVector() );
	ali1->addDiagonal( 0, 10, 20);

	for (int x = 0; x < 10; ++x)
	{
		HAlignment ali2( makeAlignmentVector() );
		ali2->addDiagonal( x, x+10, 20);
		BOOST_CHECK_EQUAL( getAlignmentIdentity( ali1, ali2, RR ), 10 - x);
		BOOST_CHECK_EQUAL( getAlignmentIdentity( ali1, ali2, CR ), 0);
		BOOST_CHECK_EQUAL( getAlignmentIdentity( ali1, ali2, RC ), 0);
		BOOST_CHECK_EQUAL( getAlignmentIdentity( ali1, ali2, CC ), 10 - x);
	}

	for (int x = 0; x < 10; ++x)
	{
		HAlignment ali2( makeAlignmentVector() );
		ali2->addDiagonal( x, x+10, 19);
		BOOST_CHECK_EQUAL( getAlignmentIdentity( ali1, ali2, RR ), 0);
		BOOST_CHECK_EQUAL( getAlignmentIdentity( ali1, ali2, CR ), 0);
		BOOST_CHECK_EQUAL( getAlignmentIdentity( ali1, ali2, RC ), 0);
		BOOST_CHECK_EQUAL( getAlignmentIdentity( ali1, ali2, CC ), 0);
	}

}

BOOST_AUTO_TEST_CASE( test_copyAlignmet1 )
{
	HAlignment ali1( makeAlignmentVector() );
	ali1->addDiagonal( 0, 10, 0);
	HAlignment dest( makeAlignmentVector() );

	for (int x = 0; x < 10; ++x)
	{
		HAlignment ali2( makeAlignmentVector() );
		ali2->addDiagonal( x, x+10, 0);

		copyAlignment( dest, ali1, ali2, RR );

		BOOST_CHECK_EQUAL(  dest->getNumAligned(), 10 - x );
	}
}


BOOST_AUTO_TEST_CASE( test_getAlignmentShortestDistance )
{
	HAlignment ali1( makeAlignmentVector() );
	ali1->addDiagonal( 20, 30, 0);
	ali1->addDiagonal( 40, 50, 0);

	{
		// adjacent before
		HAlignment ali2( makeAlignmentVector() );
		ali2->addDiagonal( 0, 10, 0);
		BOOST_CHECK_EQUAL( getAlignmentShortestDistance( ali1, ali2, RR), 11);
		BOOST_CHECK_EQUAL( getAlignmentShortestDistance( ali1, ali2, RC), 11);
		BOOST_CHECK_EQUAL( getAlignmentShortestDistance( ali1, ali2, RC), 11);
		BOOST_CHECK_EQUAL( getAlignmentShortestDistance( ali1, ali2, CC), 11);
		BOOST_CHECK_EQUAL( getAlignmentShortestDistance( ali2, ali1, RR), 11);
		BOOST_CHECK_EQUAL( getAlignmentShortestDistance( ali2, ali1, RC), 11);
		BOOST_CHECK_EQUAL( getAlignmentShortestDistance( ali2, ali1, RC), 11);
		BOOST_CHECK_EQUAL( getAlignmentShortestDistance( ali2, ali1, CC), 11);
	}

	{
		// adjacent after
		HAlignment ali2( makeAlignmentVector() );
		ali2->addDiagonal( 60, 70, 0);
		BOOST_CHECK_EQUAL( getAlignmentShortestDistance( ali1, ali2, RR), 11);
		BOOST_CHECK_EQUAL( getAlignmentShortestDistance( ali1, ali2, RC), 11);
		BOOST_CHECK_EQUAL( getAlignmentShortestDistance( ali1, ali2, RC), 11);
		BOOST_CHECK_EQUAL( getAlignmentShortestDistance( ali1, ali2, CC), 11);
		BOOST_CHECK_EQUAL( getAlignmentShortestDistance( ali2, ali1, RR), 11);
		BOOST_CHECK_EQUAL( getAlignmentShortestDistance( ali2, ali1, RC), 11);
		BOOST_CHECK_EQUAL( getAlignmentShortestDistance( ali2, ali1, RC), 11);
		BOOST_CHECK_EQUAL( getAlignmentShortestDistance( ali2, ali1, CC), 11);
	}

	{
		// inserted
		HAlignment ali2( makeAlignmentVector() );
		ali2->addDiagonal( 30, 40, 0);
		BOOST_CHECK_EQUAL( getAlignmentShortestDistance( ali1, ali2, RR), 1);
		BOOST_CHECK_EQUAL( getAlignmentShortestDistance( ali1, ali2, RC), 1);
		BOOST_CHECK_EQUAL( getAlignmentShortestDistance( ali1, ali2, RC), 1);
		BOOST_CHECK_EQUAL( getAlignmentShortestDistance( ali1, ali2, CC), 1);
		BOOST_CHECK_EQUAL( getAlignmentShortestDistance( ali2, ali1, RR), 1);
		BOOST_CHECK_EQUAL( getAlignmentShortestDistance( ali2, ali1, RC), 1);
		BOOST_CHECK_EQUAL( getAlignmentShortestDistance( ali2, ali1, RC), 1);
		BOOST_CHECK_EQUAL( getAlignmentShortestDistance( ali2, ali1, CC), 1);
	}

	{
		// overlap
		HAlignment ali2( makeAlignmentVector() );
		ali2->addDiagonal( 10, 40, 0);
		BOOST_CHECK_EQUAL( getAlignmentShortestDistance( ali1, ali2, RR), 0);
		BOOST_CHECK_EQUAL( getAlignmentShortestDistance( ali1, ali2, RC), 0);
		BOOST_CHECK_EQUAL( getAlignmentShortestDistance( ali1, ali2, RC), 0);
		BOOST_CHECK_EQUAL( getAlignmentShortestDistance( ali1, ali2, CC), 0);
		BOOST_CHECK_EQUAL( getAlignmentShortestDistance( ali2, ali1, RR), 0);
		BOOST_CHECK_EQUAL( getAlignmentShortestDistance( ali2, ali1, RC), 0);
		BOOST_CHECK_EQUAL( getAlignmentShortestDistance( ali2, ali1, RC), 0);
		BOOST_CHECK_EQUAL( getAlignmentShortestDistance( ali2, ali1, CC), 0);
	}

}



/*
BOOST_AUTO_TEST_CASE( test_fillAlignmentGaps2)
{
	HAlignment ali( makeAlignmentVector() );
	addDiagonal2Alignment( ali, 1481, 1522, -1481);
	addDiagonal2Alignment( ali, 1523, 1597, -1481);
	std::cout << AlignmentFormatEmissions( ali ) << std::endl;
	fillAlignmentGaps( ali, 3 );
	std::cout << AlignmentFormatEmissions( ali ) << std::endl;
}
*/

/*
bool TestCompressionMonotone( HAlignment & a,
		const char *xrow,
		const char *xcol,
		const char *crow = NULL,
		const char *ccol = NULL) {

	if (crow == NULL) crow = xrow;
	if (ccol == NULL) ccol = xcol;

	std::string row (xrow);
	std::string col (xcol);


	fillAlignmentCompressed( a, 3, row, 3, col);

	std::stringstream output;

	writeAlignmentCompressed( output, a  );
	output.seekp( 0);
	std::string new_row;
	std::string new_col;

	output >> new_row >> new_col;

	row = crow;
	col = ccol;

	if (row == new_row && col == new_col)
		return true;
	else {
		cout << "xrow=" << xrow << " xcol=" << xcol << endl;
		cout << "crow=" << crow << " ccol=" << ccol << endl;
		cout << *a << endl;
		cout << "row=" << row << " col=" << col << " new_row=" << new_row << " new_col=" << new_col << endl;
		return false;
	}
}

bool TestCompressionDiagonal( HAlignment & a,
		const char *xrow,
		const char *crow = NULL)
{

	if (crow == NULL) crow = xrow;

	std::string row (xrow);

	fillAlignmentCompressedDiagonal( a, row );

	std::stringstream output;
	writeAlignmentCompressedDiagonal( output, a );
	std::string new_row(output.str());

	row = crow;

	if (row == new_row)
		return true;
	else {
		cout << "xrow=" << xrow << " crow=" << crow << endl;
		cout << *a << endl;
		cout << "new_row=" << new_row << endl;
		return false;
	}
}
*/

// tests for both
void Test( HAlignment & a)
{
	{
		cout << "testing...copyAlignment()...";

		HAlignment a_new = makeAlignmentVector();
		copyAlignment( a_new, a);
	}

	{
		cout << "testing...combineAlignment()...";
		HAlignment a_new = makeAlignmentVector();
		combineAlignment( a_new, a, a, alignlib::RR);
		combineAlignment( a_new, a, a, alignlib::CR);
		combineAlignment( a_new, a, a, alignlib::RC);
		combineAlignment( a_new, a, a, alignlib::CC);
		cout << "passed" << endl;
	}

	{
		cout << "testing...rescoreAlignment()...";
		rescoreAlignment( a, 1 );
		cout << "passed" << endl;
	}
	{
		cout << "testing...rescoreAlignment()...";
		HAlignandum s1 = makeSequence("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA");
		HAlignandum s2 = makeSequence("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA");
		HScorer s( makeScorer( s1, s2 ) );
		rescoreAlignment( a, s1, s2, s);
		cout << "passed" << endl;
	}

	{
		cout << "testing...addDiagonal2Alignment()...";
		addDiagonal2Alignment( a, 1, 10, 0);
		cout << "passed" << endl;
	}

	{
		cout << "testing...addAlignment2Alignment()...";
		HAlignment b = makeAlignmentVector();
		addAlignment2Alignment( b, a );
		cout << "passed" << endl;
	}

	{
		cout << "testing...addMappedAlignment2Alignment()...";
		HAlignment b = makeAlignmentVector();
		addMappedAlignment2Alignment( b, a, a, alignlib::RR );
		cout << "passed" << endl;
	}

}

// tests for diagonal alignments
void TestDiagonal(HAlignment & a)
{

	/*
	{
		cout << "testing...fill/writeAlignmentCompressedDiagonal()...";

		if (TestCompressionDiagonal( a, "-1:-5+2" ) &&
				TestCompressionDiagonal( a, "-3:-5+2-2+2;5:-2+3" ) &&
				TestCompressionDiagonal( a, "3:-5+2" ) &&
				TestCompressionDiagonal( a, "-3:-5+2;0:-5+4" )
		)
			cout << "passed" << endl;
		else
			cout << "failed" << endl;
	}
	*/
	Test( a );
}
// tests for monotone alignments
void TestMonotone( HAlignment & a)
{
	/*
	{
		cout << "testing...fill/writeAlignmentCompressed()...";

		if (TestCompressionMonotone( a, "+5-5+5", "+15") &&				// treat gaps in one sequence correctly
				TestCompressionMonotone( a, "+15", "+5-5+5") &&				// treat gaps in one sequence correctly
				TestCompressionMonotone( a, "+5-5+5", "+5-5+5", "+10", "+10") &&	// skip a gap in both sequences
				TestCompressionMonotone( a, "-5+10", "+10-5", "+5", "+5") &&		// skip end gaps
				TestCompressionMonotone( a, "+10-5", "-5+10", "+5", "+5") &&		// skip end gaps
				TestCompressionMonotone( a, "-5+5-5+10", "-5+5-5+5-5", "+10", "+10") )	// skip end and middle gaps gaps
			cout << "passed" << endl;
		else
			cout << "failed" << endl;
	}
	 */
	{
		cout << "testing expandAlignment()...";
		a->clear();
		a->addPair( ResiduePair( 3, 4, 0));
		a->addPair( ResiduePair( 4, 5, 0));
		a->addPair( ResiduePair( 5, 7, 0));
		a->addPair( ResiduePair( 9, 9, 0));

		HAlignment b = a->getNew();
		HAlignment c = a->getNew();

		expandAlignment( b, c, a, true, true, true, true, 11, 11);
	}
	Test( a );
}
/*
int main () {

	HAlignment a;

	cout << "---------------------Testing AlignmentVector-------------------------------" << endl;
	a = makeAlignmentVector();
	TestMonotone( a );

	cout << "---------------------Testing AlignmentSet----------------------------------" << endl;
	a = makeAlignmentSet();
	TestMonotone( a );

	cout << "---------------------Testing AlignmentHash----------------------------------" << endl;
	a = makeAlignmentHash();
	TestMonotone( a );

	cout << "---------------------Testing AlignmentSetCol------------------------------" << endl;
	a = makeAlignmentSetCol();
	TestMonotone( a );

	cout << "---------------------Testing AlignmentMatrixRow-------------------------------" << endl;
	a = makeAlignmentMatrixRow();
	TestMonotone( a );

	cout << "---------------------Testing AlignmentMatrixDiagonal-------------------------------" << endl;
	a = makeAlignmentMatrixDiagonal();
	TestDiagonal( a );

	cout << "---------------------Testing AlignmentHashDiagonal------------------------------" << endl;
	a = makeAlignmentHashDiagonal();
	TestDiagonal( a );

	cout << "---------------------Testing AlignmentMatrixUnsorted-------------------------------" << endl;
	a = makeAlignmentMatrixUnsorted();
	TestDiagonal( a );

}

*/






