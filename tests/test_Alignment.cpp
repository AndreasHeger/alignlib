/*
  alignlib - a library for aligning protein sequences

  $Id: test_Alignment.cpp,v 1.4 2004/06/02 12:11:38 aheger Exp $

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
#include <cassert>

#include <time.h> 

#include "alignlib.h"
#include "Alignment.h"
#include "AlignmentIterator.h"
#include "HelpersAlignment.h"
#include "AlignlibDebug.h"

#define BOOST_TEST_MODULE
#include <boost/test/included/unit_test.hpp>
using boost::unit_test::test_suite;

using namespace std;
using namespace alignlib;

bool isIdentical( 
		const HAlignment & a, 
		const HAlignment & b, 
		bool inverse = false) 
{

	AlignmentIterator it1(a->begin());
	AlignmentIterator it1_end(a->end());

	AlignmentIterator it2(b->begin());
	AlignmentIterator it2_end(b->end());

	bool is_identical = true;

	for (; it1 != it1_end; ++it1, ++it2) 
	{
		if (!inverse) 
		{
			if (it1->mRow != it2->mRow && it1->mCol != it2->mCol) 
				is_identical = false;
		} else 
		{
			if (it1->mRow != it2->mCol && it1->mCol != it2->mRow) 
				is_identical = false;
		}
	}

	return is_identical;
}

// fill alignment with sample data
void fillAlignment( HAlignment & a)
{
	a->addPair( ResiduePair(3,3, 1.0));		
	a->addPair( ResiduePair(4,4, 1.0));				
	a->addPair( ResiduePair(5,6, 1.0));
	a->addPair( ResiduePair(6,7, 1.0));
	a->addPair( ResiduePair(8,8, 1.0));
	a->addPair( ResiduePair(9,10, 1.0));
	a->addPair( ResiduePair(10,11, 1.0));
	a->addPair( ResiduePair(12,12, 1.0));	
}
// tests for both empty and full alignments
void testAlignment( HAlignment & a)
{
	{ 
		cout << "testing...writing alignment...";
		ostringstream result;
		result << *a;
		cout << "passed" << endl;
	}

	{
		cout << "testing...iteration (post-increment)...";

		AlignmentIterator it(a->begin());
		AlignmentIterator it_end(a->end());
		for (; it != it_end; it++) 
		{ 	ResiduePair p = *it; 
			p.mRow+=1;
		}
		cout << "passed" << endl;
	}

	{
		cout << "testing...iteration (pre-increment)...";

		AlignmentIterator it(a->begin());
		AlignmentIterator it_end(a->end());
		for (; it != it_end; ++it) 
		{ 
			ResiduePair p = *it; 
			p.mRow +=1;
		}
		cout << "passed" << endl;
	}

	{ 
		cout << "testing...getNumGaps()...";
		int result = a->getNumGaps(); result++;
		cout << "passed" << endl;
	}

	{ 
		cout << "testing...getNumLength()..." ;
		int result = a->getLength(); result++;
		cout << "passed" << endl;
	}

	{ 
		cout << "testing...getScore()...";
		Score result = a->getScore(); result+=1;
		cout << "passed" << endl;
	}

	{ 
		cout << "testing...setting score...";
		a->setScore( 12.0 );
		Score score = a->getScore();
		if (score == 12.0) 
			cout << "passed" << endl;
		else
			cout << "failed" << endl;
	}

	{
		cout << "testing...getClone()...";
		HAlignment a_clone(a->getClone());

		if (a_clone->getScore()   == a->getScore() && 
				a_clone->getLength()  == a->getLength() && 
				a_clone->getNumGaps() == a->getNumGaps() &&
				isIdentical( a, a_clone ) )
			cout << "passed" << endl;
		else {
			cout << "failed" << endl;
			cout << *a << endl;
			cout << *a_clone << endl;
		}
	}

	{

		cout << "testing...getNew()...";
		HAlignment a_new(a->getNew());

		if (a_new->getScore() == 0 && 
				a_new->getLength() == 0 && 
				a_new->getNumGaps() ==0) 
			cout << "passed" << endl;
		else
			cout << "failed" << endl;

	}

	// this function assumes an ordering by row
	// which is not true for some containers.
	{ 
		cout << "testing...mapRowToCol()..."; 

		AlignmentIterator it(a->begin());
		AlignmentIterator it_end(a->end());
		Position pos = 0;
		bool passed = true;
		while (pos < a->getRowTo()) 
		{
			// check gaps
			while (it != it_end && pos < (*it).mRow)
			{
				if (a->mapRowToCol( pos++ ) != NO_POS) 
					passed = false;
			}

			if (it != it_end) 
			{
				if (a->mapRowToCol( pos) != it->mCol) 
					passed = false;
				++it;
			}
			++pos;
		}      
		if (passed)
			cout << "passed" << endl;
		else
			cout << "failed" << endl;
	}

	{ 
		cout << "testing...mapColToRow()..."; 

		AlignmentIterator it(a->begin());
		AlignmentIterator it_end(a->end());
		Position pos = 0;
		bool passed = true;
		while (pos < a->getColTo() ) 
		{
			while (it != it_end && pos < (*it).mCol) 
				if (a->mapColToRow( pos++ ) != NO_POS) passed = false;

			if (it!= it_end) 
			{
				if (a->mapColToRow( pos) != it->mRow) passed = false;
				++it;
			}

			++pos;
		}      
		if (passed)
			cout << "passed" << endl;
		else
			cout << "failed" << endl;
	}

	{ 
		cout << "testing...switchRowCol()...";

		HAlignment a_clone(a->getClone());
		a_clone->switchRowCol();    
		bool passed = isIdentical( a, a_clone, true );
		a_clone->switchRowCol();    
		passed &= isIdentical( a, a_clone );

		if (passed)
			cout << "passed" << endl;
		else
			cout << "failed" << endl;


	}

	{
		cout << "testing...removeRowRegion()..." ;

		
		if (!a->isEmpty())
		{
			HAlignment a_clone(a->getClone());

			a_clone->removeRowRegion( 3, 5);
			assert( a_clone->getRowFrom() == 5 );
			assert( a_clone->getColFrom() == 6 );
		
			a_clone->removeRowRegion( 10, 13);
			assert( a_clone->getRowTo() == 10);
			assert( a_clone->getColTo() == 11);       
			a_clone->removeRowRegion( 0, 20 );
			assert( a_clone->getLength() ==  0 );
			assert( a_clone->isEmpty()== true);
		}
        
		{
			HAlignment a_clone(a->getClone());
			a_clone->removeRowRegion(-1, -3);
			a_clone->removeRowRegion(20, 30);
			a_clone->removeRowRegion(0, 2);
			a_clone->removeRowRegion(1, 2);
			a_clone->removeRowRegion(2, 1);                               
			assert( a_clone->getLength() ==  a->getLength() );
		}
        
		cout << "passed" << endl;
	}

	{
		cout << "testing...removeColRegion()..." ;

		int i = 0;

		HAlignment  a_clone(a->getClone());

		for (i = 0; i < a->getColTo() + 5; i+=3) 
			a_clone->removeColRegion(i, i+3);       

		cout << "passed" << endl;
	}

	// test removing all pairs
	{
		HAlignment a_clone(a->getClone());

		AlignmentIterator it(a->begin());
		AlignmentIterator it_end(a->end());
		
		int naligned = a_clone->getNumAligned();
		for (; it != it_end; ++it)
		{
			a_clone->removePair( *it );
			BOOST_CHECK_EQUAL( a_clone->getNumAligned(), --naligned );
		}
		BOOST_CHECK_EQUAL(a_clone->getLength(), 0 );
		BOOST_CHECK_EQUAL(a_clone->getNumGaps() , 0);
		BOOST_CHECK_EQUAL(a_clone->getNumAligned() , 0);
		BOOST_CHECK_EQUAL(a_clone->isEmpty(), true );
	}
	
	{ 
		HAlignment a_clone( a->getClone());
		a_clone->clear();
		BOOST_CHECK_EQUAL(a_clone->getScore(), 0);
		BOOST_CHECK_EQUAL(a_clone->getLength(), 0 );
		BOOST_CHECK_EQUAL(a_clone->getNumGaps() , 0);
		BOOST_CHECK_EQUAL(a_clone->getNumAligned() , 0);		
		BOOST_CHECK_EQUAL(a_clone->isEmpty(), true );
	}

	if ( a->isEmpty() )
		return;
	
	BOOST_CHECK_EQUAL( (*(a->begin())).mRow, a->getRowFrom() ); 	
	BOOST_CHECK_EQUAL( a->front().mRow, a->getRowFrom() ); 		

	BOOST_CHECK_EQUAL( (*(a->begin())).mCol, a->getColFrom() ); 	
	BOOST_CHECK_EQUAL( a->front().mCol, a->getColFrom() );


}

//----------------------------------------------------------
// main test routine for a pairwise alignment
void runTests( HAlignment & a ) 
{
	testAlignment( a );
	fillAlignment( a );
	testAlignment( a );
	a->clear();
	fillAlignment( a );
	testAlignment( a );
}

#define create_test( name, factory ) \
	BOOST_AUTO_TEST_CASE( name ) { HAlignment a(factory()); runTests(a); }

create_test( Vector, makeAlignmentVector );
create_test( Set, makeAlignmentSet );
create_test( Hash, makeAlignmentHash );
create_test( HashDiagonal, makeAlignmentHashDiagonal );
create_test( SetCol, makeAlignmentSetCol );
create_test( MatrixRow, makeAlignmentMatrixRow );
create_test( MatrixDiagonal, makeAlignmentMatrixDiagonal );
create_test( MatrixUnsorted, makeAlignmentMatrixUnsorted );



