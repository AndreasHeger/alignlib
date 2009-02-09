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

#include "Alignandum.h"
#include "HelpersAlignandum.h"
#include "Alignator.h"
#include "HelpersAlignator.h"
#include "Alignatum.h"
#include "HelpersAlignatum.h"
#include "Alignment.h"
#include "HelpersAlignment.h"
#include "MultAlignment.h"
#include "HelpersMultAlignment.h"
#include "AlignlibDebug.h"
#include "AlignlibException.h"

using namespace std;
using namespace alignlib;

const char FILE_SEQ[] = "data/test.seq";
const char FILE_FASTA[] = "data/test.fasta";

#define BOOST_TEST_MAIN
#include <boost/test/included/unit_test.hpp>
using boost::unit_test::test_suite;

class my_exception{};

void test_GenericMultAlignment(
		const HMultAlignment & r )
{
	int nseqs = 10;

	// checking an empty alignment
	{
		HMultAlignment clone(r->getNew());
		BOOST_CHECK_EQUAL( clone->getLength(), 0);
		BOOST_CHECK_EQUAL( clone->getNumSequences(), 0);
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

}

BOOST_AUTO_TEST_CASE( test_MultAlignment )
{
	HMultAlignment r( makeMultAlignment());
	test_GenericMultAlignment( r );
}




