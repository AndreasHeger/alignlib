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
#include "alignlib_fwd.h"

#include "Alignandum.h"
#include "Encoder.h"
#include "HelpersEncoder.h"
#include "HelpersAlignandum.h"

using namespace std;

using namespace alignlib;

void checkingStart( const std::string & s)
{
	// std::cout << "starting check:" << s << "...";
}

void checkingEnd( bool passed = true)
{
	if (passed)
	{}
		// std::cout << "passed" << std::endl;
	else
	{}
		//std::cout << "failed" << std::endl;
}

void runTests( HAlignandum & a, const std::string & sample )
{

	const HEncoder translator = a->getToolkit()->getEncoder();

	checkingStart( "translation and range" );
	{
	assert( a->getFrom() == 0);
	assert( a->getTo() == sample.size() );
	assert( sample.size() == a->getLength() );
	for ( int x = 0; x < sample.size(); ++x)
		assert( a->asChar(x) == sample[x]);

	assert( a->asString() == sample );
	}
	checkingEnd();

	// check that useSegment does not interfere with output
	checkingStart( "segments" );
	{
		HAlignandum clone(a->getClone());
		Position from = 0;
		Position to = sample.size();
		clone->useSegment( from, to);
		assert( clone->getFrom() == from);
		assert( clone->getTo() == to );
		for ( int x = 0; x < sample.size(); ++x)
		{
			assert( clone->asChar(x) == sample[x]);
		}
		assert( clone->asString() == sample );
	}
	checkingEnd();

	checkingStart( "saving" );
	// check saving/loading from stream.
	{
		ofstream file("test_Alignandum.tmp", ios::binary);
		a->save( file );
		a->save( file );
		file.close();
	}
	checkingEnd();

	checkingStart( "loading" );
	{
		ifstream file("test_Alignandum.tmp", ios::binary) ;

		HAlignandum b;
		int n = 0;
		while ( file.peek() != EOF )
		{
			b = loadAlignandum( file );
			assert( a->getFrom() == b->getFrom() );
			assert( a->getLength() == b->getLength() );
			assert( a->getTo() == a->getTo() );
			assert( a->asString() == b->asString() );
			++n;
		}
		assert( n == 2 );
	}
	checkingEnd();

	// check masking
	checkingStart( "masking" );
	{
		HAlignandum clone(a->getClone());

		clone->mask( 0, a->getLength() );
		for (Position p = 0; p < clone->getLength(); ++p)
		{
			assert( clone->asResidue(p) == translator->getMaskCode() );
			assert( clone->asChar(p) == translator->getMaskChar() );
		}
	}
	checkingEnd();

	checkingStart( "shuffling" );
	{
		HAlignandum clone(a->getClone());
		clone->shuffle();
		assert( clone->getLength() == a->getLength() );
	}
	checkingEnd();

	{
		HAlignandum clone(a->getClone());
		a->setStorageType( Full );
		checkingStart( "writing" );
		{
			std::ofstream outfile;
			outfile.open( "test_full.out", std::ios::binary );
			a->save( outfile );
			outfile.close();
		}
		checkingEnd();

		checkingStart( "reading" );
		{
			std::ifstream infile;
			infile.open( "test_full.out", std::ios::binary );
			HAlignandum a(loadAlignandum( infile ) );
			infile.close();
		}
		checkingEnd();
	}

	{
		HAlignandum clone(a->getClone());
		a->setStorageType( Sparse );
		checkingStart( "writing" );
		{
			std::ofstream outfile;
			outfile.open( "test_sparse.out", std::ios::binary );
			a->save( outfile );
			outfile.close();
		}
		checkingEnd();

		checkingStart( "reading" );
		{
			std::ifstream infile;
			infile.open( "test_sparse.out", std::ios::binary );
			HAlignandum a(loadAlignandum( infile ) );
			infile.close();
		}
		checkingEnd();
	}

}

void testAlignandum( HAlignandum & a, const std::string & sample )
{
	// std::cout << "--- testing fresh --- " << std::endl;
	runTests( a, sample);
	// std::cout << "--- testing prepared --- " << std::endl;
	a->prepare();
	runTests( a, sample);
	// std::cout << "--- testing released --- " << std::endl;
	a->release();
	runTests( a, sample);
}


int main ()
{

	std::string ref_protein20 = "ACDEFGHIKLMNPQRSTVWY";
	std::string ref_protein20x3 = ref_protein20 + ref_protein20 + ref_protein20;
	{
		// std::cout << "--- testing Sequence ----" << std::endl;
		HAlignandum a(makeSequence( "ACA") );
		testAlignandum( a, "ACA" );
	}


	{
		//std::cout << "--- testing masked Sequence ----" << std::endl;
		HAlignandum a(makeSequence( "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX") );
		testAlignandum( a, "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" );
	}

	{
		// std::cout << "--- testing Sequence - protein 20 ----" << std::endl;
		HAlignandum a(makeSequence( ref_protein20 ) );
		testAlignandum( a, ref_protein20 );
	}

	{
		// std::cout << "--- testing Profile----" << std::endl;
		// HAlignandum a(makeProfile( ref_protein20x3, 3) );
		// testAlignandum( a, ref_protein20  );
		HAlignandum a(makeProfile( std::string("AAA"), 1) );
		testAlignandum( a, std::string("AAA")  );
	}

	return EXIT_SUCCESS;

}








