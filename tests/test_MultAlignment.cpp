/*
  alignlib - a library for aligning protein sequences

  $Id: test_MultAlignment.cpp,v 1.6 2004/06/02 12:11:38 aheger Exp $

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

/** Test the MultAlignment - object
 */

#if HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>
#include <fstream>
#include <sstream>
#include <cassert>

#include "Matrix.h"
#include "Alignandum.h"
#include "HelpersAlignandum.h"
#include "Alignator.h"
#include "HelpersAlignator.h"
#include "Alignatum.h"
#include "HelpersAlignatum.h"
#include "Alignment.h"
#include "HelpersAlignment.h"
#include "MultAlignment.h"
#include "MultAlignmentFormat.h"
#include "HelpersMultAlignment.h"
#include "AlignlibDebug.h"
#include "AlignlibException.h"
#include "MultipleAlignment.h"
#include "HelpersMultipleAlignment.h"

using namespace std;
using namespace alignlib;

const char FILE_SEQ[] = "data/test.seq";
const char FILE_FASTA[] = "data/test.fasta";

#define BOOST_TEST_MAIN
#include <boost/test/included/unit_test.hpp>
using boost::unit_test::test_suite;

class my_exception{};

// simple test of expansion: only one sequence without gaps
void test_Expand0( const HMultAlignment & r )
{
	{
		HMultAlignment clone(r->getNew());
		{
			HAlignment ali(makeAlignmentVector());
			ali->addDiagonal( 0,8,0);
			clone->add( ali );
		}
		HMultAlignment clone2(clone->getClone());
		clone2->expand( HAlignandumVector(new AlignandumVector()));
		BOOST_CHECK_EQUAL( checkMultAlignmentIdentity( clone, clone2), true);
	}

	// check for double expansion
	{
		HMultAlignment clone(r->getNew());
		{
			HAlignment ali(makeAlignmentVector());
			ali->addDiagonal( 0,4,0);
			ali->addDiagonal( 4,8,4);
			clone->add( ali );
		}
		HMultAlignment clone2(clone->getClone());
		clone2->expand( HAlignandumVector(new AlignandumVector()));
		BOOST_CHECK_EQUAL( checkMultAlignmentIdentity( clone, clone2), false);
		clone->expand( HAlignandumVector(new AlignandumVector()));
		BOOST_CHECK_EQUAL( checkMultAlignmentIdentity( clone, clone2), true);
		clone->expand( HAlignandumVector(new AlignandumVector()));
		BOOST_CHECK_EQUAL( checkMultAlignmentIdentity( clone, clone2), true);
	}
}

// test expansion of multiple alignment
void test_Expand1(
		const HMultAlignment & r )
{

	HMultAlignment copy(r->getNew());
	{

		{
			HAlignment ali(makeAlignmentVector());
			ali->addDiagonal( 0,3,+2 );
			ali->addDiagonal( 3,6,+4 );
			copy->add( ali );
		}
		{
			HAlignment ali(makeAlignmentVector());
			ali->addDiagonal( 0,1,+1 );
			ali->addDiagonal( 1,6,+3 );
			copy->add( ali );
			copy->add( ali );
		}
	}

	// Input:
	// 012345
	// 234789
	// 145678
	// 145678

	// expansion within mali
	// Output:
	// 0123456789 10 11
	// 2----34567 8  9
	// 123--45--6 7  8
	// 1--2345--6 7  8
	{
		HMultAlignment clone(copy->getClone());

		clone->expand( HAlignandumVector(new AlignandumVector()));
		{
			HAlignment ali_test(makeAlignmentVector());
			ali_test->addDiagonal( 0,1,+2 );
			ali_test->addDiagonal( 5,12,-2 );
			BOOST_CHECK_EQUAL( checkAlignmentIdentity( ali_test, clone->getRow(0)), true);
		}

		{
			HAlignment ali_test(makeAlignmentVector());
			ali_test->addDiagonal( 0,3,+1 );
			ali_test->addDiagonal( 5,7,-1 );
			ali_test->addDiagonal( 9,12,-3 );
			BOOST_CHECK_EQUAL( checkAlignmentIdentity( ali_test, clone->getRow(1)), true);
		}
		{
			HAlignment ali_test(makeAlignmentVector());
			ali_test->addDiagonal( 0,1,+1 );
			ali_test->addDiagonal( 3,7,-1 );
			ali_test->addDiagonal( 9,12,-3 );
			BOOST_CHECK_EQUAL( checkAlignmentIdentity( ali_test, clone->getRow(2)), true);
		}
		{
			HMultAlignment clone2(clone->getClone());
			clone2->expand( HAlignandumVector(new AlignandumVector()));
			BOOST_CHECK_EQUAL( checkMultAlignmentIdentity( clone, clone2), true);
		}

		// check if collapse gets us back to start
		{
			clone->shrink();
			BOOST_CHECK_EQUAL( checkMultAlignmentIdentity( clone, copy), true);
		}

		// check if collapse is not over-eager
		{
			clone->shrink();
			BOOST_CHECK_EQUAL( checkMultAlignmentIdentity( clone, copy), true);
		}


	}

	// expansion within and outside of mali
	// Output:
	// 0123456789 10 11 12 13 14 15 16 17
  	// 01--2----3 4  5   6  7  8  9  -  -
	// --0-123--4 5  -   -  6  7  8  9  -
	// ---01--234 5  -   -  6  7  8  -  9
	{
		HMultAlignment clone(copy->getClone());

		HAlignandumVector seqs(new AlignandumVector());
		seqs->push_back( makeSequence( "AAAAAAAAAA" ) );
		seqs->push_back( makeSequence( "AAAAAAAAAA" ) );
		seqs->push_back( makeSequence( "AAAAAAAAAA" ) );

		clone->expand( seqs );
		{
			HAlignment ali_test(makeAlignmentVector());
			ali_test->addDiagonal( 0,2, 0 );
			ali_test->addDiagonal( 4,5,-2 );
			ali_test->addDiagonal( 9,16,-6 );
			BOOST_CHECK_EQUAL( checkAlignmentIdentity( ali_test, clone->getRow(0)), true);
		}

		{
			HAlignment ali_test(makeAlignmentVector());
			ali_test->addDiagonal( 2, 3, -2 );
			ali_test->addDiagonal( 4, 7, -3 );
			ali_test->addDiagonal( 9,11, -5 );
			ali_test->addDiagonal(13,17, -7 );
			BOOST_CHECK_EQUAL( checkAlignmentIdentity( ali_test, clone->getRow(1)), true);
		}
		{
			HAlignment ali_test(makeAlignmentVector());
			ali_test->addDiagonal( 3,5,-3 );
			ali_test->addDiagonal( 7,11,-5 );
			ali_test->addDiagonal( 13,16,-7 );
			ali_test->addDiagonal( 17,18,-8 );
			BOOST_CHECK_EQUAL( checkAlignmentIdentity( ali_test, clone->getRow(2)), true);
		}
		{
			HMultAlignment clone2(clone->getClone());
			clone2->expand( seqs );
			BOOST_CHECK_EQUAL( checkMultAlignmentIdentity( clone, clone2), true);
		}
		// check if collapse gets us back to start
		{
			clone->shrink();
			BOOST_CHECK_EQUAL( checkMultAlignmentIdentity( clone, copy), true);
		}

		// check if collapse is not over-eager
		{
			clone->shrink();
			BOOST_CHECK_EQUAL( checkMultAlignmentIdentity( clone, copy), true);
		}
	}
}

// test for expansion: unaligned sequences
void test_Expand2( const HMultAlignment & r )
{
	HMultAlignment clone(r->getNew());
	{
		HAlignment ali(makeAlignmentVector());
		ali->addDiagonal( 0,8,0);
		clone->add( ali );
		clone->add( makeAlignmentVector() );
	}
	// expansion without sequences: no change to mali
	{
		HMultAlignment clone2(clone->getClone());
		clone2->expand( HAlignandumVector(new AlignandumVector()));
		BOOST_CHECK_EQUAL( checkMultAlignmentIdentity( clone, clone2), true);
	}

	// expansion with sequences
	{
		HAlignandumVector seqs(new AlignandumVector());
		seqs->push_back( makeSequence( "AAAAAAAA" ) );
		seqs->push_back( makeSequence( "AAAAAAAA" ) );
		HMultAlignment clone2(clone->getClone());
		clone2->expand( seqs );
		BOOST_CHECK_EQUAL( checkAlignmentIdentity( (*clone)[0], (*clone2)[0]), true);
		BOOST_CHECK_EQUAL( (*clone2)[0]->getColFrom(), (*clone2)[1]->getColFrom());
		BOOST_CHECK_EQUAL( (*clone2)[0]->getColTo(), (*clone2)[1]->getColTo());
		BOOST_CHECK_EQUAL( 2 * (*clone2)[0]->getRowTo(), (*clone2)[1]->getRowTo());
		BOOST_CHECK_EQUAL( (*clone2)[0]->getRowTo(), (*clone2)[1]->getRowFrom());
	}
}

// test for expansion: empty alignment
void test_Expand3( const HMultAlignment & r )
{
	HMultAlignment clone(r->getNew());
	{
		clone->add( makeAlignmentVector() );
		clone->add( makeAlignmentVector() );
		BOOST_CHECK_EQUAL( clone->isEmpty(), false);
		BOOST_CHECK_EQUAL( clone->getNumSequences(), 2);
		BOOST_CHECK_EQUAL( clone->getLength(), 0);
	}

	// expansion without sequences: no change to mali
	{
		HMultAlignment clone2(clone->getClone());
		clone2->expand( HAlignandumVector(new AlignandumVector()));
		BOOST_CHECK_EQUAL( clone2->isEmpty(), false);
		BOOST_CHECK_EQUAL( clone2->getNumSequences(), 2);
		BOOST_CHECK_EQUAL( clone2->getLength(), 0);
	}

	// expansion with sequences
	{
		HAlignandumVector seqs(new AlignandumVector());
		seqs->push_back( makeSequence( "AAAAAAAA" ) );
		seqs->push_back( makeSequence( "AAAAAAAA" ) );
		HMultAlignment clone2(clone->getClone());
		clone2->expand( seqs );
		BOOST_CHECK_EQUAL( clone2->isEmpty(), false);
		BOOST_CHECK_EQUAL( clone2->getNumSequences(), 2);
		BOOST_CHECK_EQUAL( clone2->getLength(), 16);
		BOOST_CHECK_EQUAL( (*clone2)[0]->getRowFrom(), 0);
		BOOST_CHECK_EQUAL( (*clone2)[0]->getRowTo(), 8);
		BOOST_CHECK_EQUAL( (*clone2)[0]->getColFrom(), 0);
		BOOST_CHECK_EQUAL( (*clone2)[0]->getColTo(), 8);

	}
}

// test if alignment works with wrap-around alignments
void test_WrapAround( const HMultAlignment & r )
{
	return;
	HMultAlignment clone(r->getNew());
	{
		HAlignment ali = makeAlignmentVector();
		ali->addDiagonal( 0, 10, 10 );
		ali->addDiagonal( 10, 20, -10 );
		std::cout << *ali << std::endl;
		clone->add( ali );
		BOOST_CHECK_EQUAL( clone->isEmpty(), false);
		BOOST_CHECK_EQUAL( clone->getNumSequences(), 1);
		BOOST_CHECK_EQUAL( clone->getLength(), 10);
		std::cout << *clone << std::endl;
	}

}

void fillMali( HMultAlignment & mali, const HAlignment & a )
{
	mali->clear();
	{
		HAlignment ali(a->getNew());
		ali->addDiagonal( 0,3,+2 );
		ali->addDiagonal( 3,6,+4 );
		mali->add( ali );
	}
	{
		HAlignment ali(a->getNew());
		ali->addDiagonal( 0,1,+1 );
		ali->addDiagonal( 1,6,+3 );
		mali->add( ali );
		mali->add( ali );
	}
}

void test_GenericMultAlignment(
		const HMultAlignment & r )
{
	int nseqs = 10;

	// checking an empty alignment
	{
		HMultAlignment clone(r->getNew());
		BOOST_CHECK_EQUAL( clone->getLength(), 0);
		BOOST_CHECK_EQUAL( clone->getNumSequences(), 0);
		BOOST_CHECK_EQUAL( clone->getLength(), clone->getColumnCounts()->size() );
		BOOST_CHECK_EQUAL( clone->getNumSequences(), clone->getRowCounts()->size() );
	}

	// checking out-of-range access of an empty alignment
	{
		HMultAlignment clone(r->getNew());
		BOOST_CHECK_THROW( (*clone)[0],  AlignlibException);
		BOOST_CHECK_THROW( (*clone)[1],  AlignlibException);
		BOOST_CHECK_THROW( clone->getRow(0),  AlignlibException);
		BOOST_CHECK_THROW( clone->getRow(1),  AlignlibException);
	}

	// check adding without gaps and out-of-range access
	{
		HMultAlignment clone(r->getNew());
		for (int x = 0; x < nseqs; ++x)
		{
			HAlignment ali(makeAlignmentVector());
			ali->addDiagonal( 0, nseqs, 0);
			clone->add(ali);
			BOOST_CHECK_EQUAL( checkAlignmentIdentity( (*clone)[x], ali ), true );
		}
		BOOST_CHECK_EQUAL( clone->getLength(), nseqs);
		BOOST_CHECK_EQUAL( clone->getNumSequences(), nseqs);
		BOOST_CHECK_THROW( (*clone)[nseqs+1], AlignlibException);
		BOOST_CHECK_THROW( (*clone)[-1], AlignlibException);
		for (int x = 0; x < nseqs; ++x)
			BOOST_CHECK_EQUAL( clone->isAligned(x), true );

		BOOST_CHECK_EQUAL( clone->getLength(), clone->getColumnCounts()->size() );
		BOOST_CHECK_EQUAL( clone->getNumSequences(), clone->getRowCounts()->size() );
		{
			HCountVector counts(clone->getColumnCounts());
			for (int x = 0; x < counts->size(); ++x)
				BOOST_CHECK_EQUAL( (*counts)[x], nseqs );
		}
		{
			HCountVector counts(clone->getRowCounts());
			for (int x = 0; x < counts->size(); ++x)
				BOOST_CHECK_EQUAL( (*counts)[x], nseqs );
		}

	}

	// check adding with gaps, out-of-range access and aligned columns
	{
		HMultAlignment clone(r->getNew());
		for (int x = 0; x < nseqs; ++x)
		{
			HAlignment ali(makeAlignmentVector());
			ali->addDiagonal( 0, nseqs + x, 0);
			ali->addDiagonal( 2 * nseqs, 3 * nseqs, 0);
			clone->add(ali);
			BOOST_CHECK_EQUAL( checkAlignmentIdentity( (*clone)[x], ali ), true );
			for (int y = 0; y < nseqs; ++y) BOOST_CHECK_EQUAL( clone->isAligned(y), true );
			for (int y = nseqs; y < nseqs + x; ++y) BOOST_CHECK_EQUAL( clone->isAligned(y), true );
			for (int y = nseqs + x; y < 2 * nseqs; ++y) BOOST_CHECK_EQUAL( clone->isAligned(y), false );
			for (int y = 2 * nseqs; y < 3 * nseqs; ++y) BOOST_CHECK_EQUAL( clone->isAligned(y), true );

		}
		BOOST_CHECK_EQUAL( clone->getLength(), 3 * nseqs);
		BOOST_CHECK_EQUAL( clone->getNumSequences(), nseqs);
		BOOST_CHECK_THROW( (*clone)[nseqs+1], AlignlibException);
		BOOST_CHECK_THROW( (*clone)[-1], AlignlibException);

		BOOST_CHECK_EQUAL( clone->getLength(), clone->getColumnCounts()->size() );
		BOOST_CHECK_EQUAL( clone->getNumSequences(), clone->getRowCounts()->size() );
		{
			HCountVector counts(clone->getColumnCounts());
			for (int x = 0; x < nseqs; ++x)
				BOOST_CHECK_EQUAL( (*counts)[x], nseqs );
			for (int x = 0; x < nseqs; ++x)
				BOOST_CHECK_EQUAL( (*counts)[nseqs+x], nseqs - x - 1);
			for (int x = 2*nseqs; x < 3*nseqs; ++x)
				BOOST_CHECK_EQUAL( (*counts)[x], nseqs );
		}
		{
			HCountVector counts(clone->getRowCounts());
			for (int x = 0; x < counts->size(); ++x)
				BOOST_CHECK_EQUAL( (*counts)[x], nseqs + nseqs + x);
		}

	}

	// check adding mali to mali with gaps, out-of-range access and aligned columns
	{
		HMultAlignment clone1(r->getNew());
		HMultAlignment clone2(r->getNew());
		// build alignment of the form
		// xxxxxx
		// xxxxxx
		for (int x = 0; x < nseqs; ++x)
		{
			HAlignment ali(makeAlignmentVector());
			ali->addDiagonal( 0, nseqs, 0);
			clone1->add(ali);
		}
		// build alignment of the form
		// xxxxx---
		// -xxxxx--
		// --xxxxx-
		// ---xxxxx
		for (int x = 0; x < nseqs; ++x)
		{
			HAlignment ali(makeAlignmentVector());
			ali->addDiagonal( x, nseqs + x, 0);
			clone2->add(ali);
		}

		HAlignment map_b2a(makeAlignmentVector());
		int half = (int)(nseqs / 2);
		map_b2a->addDiagonal( 0, half, 0);
		map_b2a->addDiagonal( nseqs, nseqs + half, -half );

		clone2->add( clone1, map_b2a );
		BOOST_CHECK_EQUAL( clone2->getLength(), 2 * nseqs - 1);
		BOOST_CHECK_EQUAL( clone2->getNumSequences(), 2 * nseqs);
		BOOST_CHECK_EQUAL( checkAlignmentIdentity( map_b2a, (*clone2)[nseqs+1] ), true );
	}

	// check adding mali to mali with gaps, out-of-range access and aligned columns
	{
		HMultAlignment clone1(r->getNew());
		HMultAlignment clone2(r->getNew());
		// build alignment of the form
		// xxxxxxxxxx
		// xxxxxxxxxx
		for (int x = 0; x < nseqs; ++x)
		{
			HAlignment ali(makeAlignmentVector());
			ali->addDiagonal( 0, nseqs, 0);
			clone1->add(ali);
		}
		// build alignment of the form
		// xxxxxxxxxx---
		// -xxxxxxxxxx--
		// --xxxxxxxxxx-
		// ---xxxxxxxxxx
		for (int x = 0; x < nseqs; ++x)
		{
			HAlignment ali(makeAlignmentVector());
			ali->addDiagonal( x, nseqs + x, 0);
			clone2->add(ali);
		}

		HAlignment map_old2new(makeAlignmentVector());
		int half = (int)(nseqs / 2);
		map_old2new->addDiagonal( 0, half, 0);
		map_old2new->addDiagonal( nseqs, nseqs + half, -half );

		clone2->add( clone1, map_old2new, map_old2new );
		BOOST_CHECK_EQUAL( clone2->getLength(), nseqs);
		BOOST_CHECK_EQUAL( clone2->getNumSequences(), 2 * nseqs);
		for (int x = 0; x < map_old2new->getColTo(); ++x)
		{
			Position pos = map_old2new->mapColToRow( x );
			if (pos == NO_POS)
				BOOST_CHECK_EQUAL( clone2->isAligned( x ), false );
			else
			{
				BOOST_CHECK_EQUAL( clone2->isAligned( x ), true );
			}
		}
		HPositionMatrix matrix(clone2->getPositionMatrix());
	}

	test_Expand0( r );
	test_Expand1( r );
	test_Expand2( r );
	test_Expand3( r );

	// test wrap-around alignments
	test_WrapAround( r );

}

BOOST_AUTO_TEST_CASE( test_MultAlignment_Blocks )
{
	HMultAlignment r( makeMultAlignment() );
	fillMali( r, makeAlignmentBlocks() );
	r->getGapCounts( HAlignandumVector(new AlignandumVector() ), AggCount );
}


BOOST_AUTO_TEST_CASE( test_MultAlignment )
{
	HMultAlignment r( makeMultAlignment());
	test_GenericMultAlignment( r );
}

BOOST_AUTO_TEST_CASE( test_Conversion )
{
	HMultipleAlignment reference(makeMultipleAlignment());

	{
		reference->add(makeAlignatum("0123456789"));
		reference->add(makeAlignatum("0123456789"));
		reference->add(makeAlignatum("0123456789"));
	}

	HMultAlignment r( makeMultAlignment() );
	HAlignment a(makeAlignmentVector());
	a->addDiagonal( 0, reference->getLength(), 0);
	r->add( reference, a );
	BOOST_CHECK_EQUAL( r->getLength(), reference->getLength() );
	BOOST_CHECK_EQUAL( r->getNumSequences(), reference->getNumSequences() );

	for (Position x = 0; x < 10; ++x)
	{
		BOOST_CHECK_EQUAL( r->isAligned( x ), true );
	}

}

BOOST_AUTO_TEST_CASE( test_ConversionWithMapping )
{
	HMultipleAlignment reference(makeMultipleAlignment());

	{
		reference->add(makeAlignatum("0123456789"));
		reference->add(makeAlignatum("0123456789"));
		reference->add(makeAlignatum("0123456789"));
	}

	HMultAlignment r( makeMultAlignment() );
	HAlignment a(makeAlignmentVector());
	a->addDiagonal( 0, 5, 0);
	a->addDiagonal( 5, 10, 1);
	r->add( reference, a );
	BOOST_CHECK_EQUAL( r->getLength()-1, reference->getLength() );
	BOOST_CHECK_EQUAL( r->getNumSequences(), reference->getNumSequences() );

	for (Position x = 0; x < 5; ++x)
		BOOST_CHECK_EQUAL( r->isAligned( x ), true );
	BOOST_CHECK_EQUAL( r->isAligned( 5 ), false );
	for (Position x = 6; x < 11; ++x)
		BOOST_CHECK_EQUAL( r->isAligned( x ), true );
}

BOOST_AUTO_TEST_CASE( test_ConversionGap )
{
	HMultipleAlignment reference(makeMultipleAlignment());

	{
		reference->add(makeAlignatum("01234-6789"));
		reference->add(makeAlignatum("01234-6789"));
		reference->add(makeAlignatum("01234-6789"));
	}

	HMultAlignment r( makeMultAlignment() );
	HAlignment a(makeAlignmentVector());
	a->addDiagonal( 0, reference->getLength(), 0);
	r->add( reference, a );
	BOOST_CHECK_EQUAL( r->getLength(), reference->getLength() );
	BOOST_CHECK_EQUAL( r->getNumSequences(), reference->getNumSequences() );

	for (Position x = 0; x < 5; ++x)
		BOOST_CHECK_EQUAL( r->isAligned( x ), true );
	BOOST_CHECK_EQUAL( r->isAligned( 5 ), false );
	for (Position x = 6; x < 10; ++x)
		BOOST_CHECK_EQUAL( r->isAligned( x ), true );

}

BOOST_AUTO_TEST_CASE( test_ConversionGaps )
{
	HMultipleAlignment reference(makeMultipleAlignment());

	{
		reference->add(makeAlignatum("-1234-6789-"));
		reference->add(makeAlignatum("-1234-6789-"));
		reference->add(makeAlignatum("-1234-6789-"));
	}

	HMultAlignment r( makeMultAlignment() );
	HAlignment a(makeAlignmentVector());
	a->addDiagonal( 0, reference->getLength(), 0);
	r->add( reference, a );
	BOOST_CHECK_EQUAL( r->getLength(), reference->getLength() );
	BOOST_CHECK_EQUAL( r->getNumSequences(), reference->getNumSequences() );

	BOOST_CHECK_EQUAL( r->isAligned( 0 ), false );
	for (Position x = 1; x < 5; ++x)
		BOOST_CHECK_EQUAL( r->isAligned( x ), true );
	BOOST_CHECK_EQUAL( r->isAligned( 5 ), false );
	for (Position x = 6; x < 10; ++x)
		BOOST_CHECK_EQUAL( r->isAligned( x ), true );
	BOOST_CHECK_EQUAL( r->isAligned( 10 ), false );
}

