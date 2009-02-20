/*
  alignlib - a library for aligning protein sequences

  $Id: test_Alignator.cpp,v 1.3 2004/06/02 12:11:38 aheger Exp $

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

    note: for AlignatorIdentity, AlignatorSimilarity, etc. the
    alignments look terrible, since they are wraparound pairwise
    alignments.

 */

#if HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

#include <time.h>

#include "alignlib.h"
#include "MultipleAlignator.h"
#include "MultAlignment.h"
#include "Alignator.h"
#include "HelpersAlignator.h"

#define BOOST_TEST_MODULE
#include <boost/test/included/unit_test.hpp>
using boost::unit_test::test_suite;

using namespace std;
using namespace alignlib;

void testMultipleAlignator(
		HMultipleAlignator & ma,
		std::vector< std::string > & result1,
		std::vector< std::string > & result2,
		bool first_is_empty )
{
	HMultAlignment result(makeMultAlignment());

	{
		HStringVector sequences(new StringVector());

		// test empty alignment
		{
			ma->align( result, sequences );
			BOOST_CHECK_EQUAL( result->getNumSequences(), 0);
			BOOST_CHECK_EQUAL( result->getLength(), 0);
		}

		// add some sequences
		sequences->push_back( "EEEEAAAADDDDMMMMEEEE");
		sequences->push_back( "FFFFAAAACCCCMMMMFFFF");
		sequences->push_back( "GGGGAAAAKKKKMMMMGGGG");
		{
			ma->align( result, sequences );
			BOOST_CHECK_EQUAL( result->getNumSequences(), sequences->size());
			MultAlignmentFormatPlain f( result, sequences);
			for (int x = 0; x < sequences->size(); ++x)
				BOOST_CHECK_EQUAL( f.mData[x]->getString(), result1[x] );
		}

		// test adding empty sequence
		sequences->push_back( "" );
		{
			ma->align( result, sequences );
			BOOST_CHECK_EQUAL( result->getNumSequences(), sequences->size());
			MultAlignmentFormatPlain f( result, sequences);
			for (int x = 0; x < sequences->size() - 1; ++x)
				BOOST_CHECK_EQUAL( f.mData[x]->getString(), result1[x] );
			BOOST_CHECK_EQUAL( f.mData[sequences->size()-1]->getString(), "");
		}

	// test adding empty sequence
		sequences->insert( sequences->begin(), "" );
		{
			ma->align( result, sequences );
			BOOST_CHECK_EQUAL( result->getNumSequences(), sequences->size());
			MultAlignmentFormatPlain f( result, sequences);
			if (!first_is_empty)
			{
				BOOST_CHECK_EQUAL( f.mData[0]->getString(), "");
				for (int x = 1; x < sequences->size() - 1; ++x)
					BOOST_CHECK_EQUAL( f.mData[x]->getString(), result1[x-1] );
				BOOST_CHECK_EQUAL( f.mData[sequences->size()-1]->getString(), "");
			}
			else
			{
				BOOST_CHECK_EQUAL( result->getLength(), 0);
				for (int x = 0; x < sequences->size(); ++x)
					BOOST_CHECK_EQUAL( f.mData[x]->getString(), "");
			}
		}
	}

	// test using ranges
	{
		HAlignandumVector sequences(new AlignandumVector());
		sequences->push_back( makeSequence("EEEEAAAADDDDMMMMEEEE") );
		sequences->push_back( makeSequence("FFFFAAAACCCCMMMMFFFF") );
		sequences->push_back( makeSequence("GGGGAAAAKKKKMMMMGGGG") );
		for (int x = 0; x < sequences->size(); ++x)
			(*sequences)[x]->useSegment( 2, 18 );

		{
			ma->align( result, sequences );
			BOOST_CHECK_EQUAL( result->getNumSequences(), sequences->size());
			MultAlignmentFormatPlain f( result, sequences);
			for (int x = 0; x < sequences->size(); ++x)
				BOOST_CHECK_EQUAL( f.mData[x]->getString(), result2[x] );
		}

	}

}

BOOST_AUTO_TEST_CASE( multiple_alignment_simple )
{
	setDefaultSubstitutionMatrix( makeSubstitutionMatrix( getDefaultEncoder()->getAlphabetSize(), 1, -1000) );
	HAlignator a(makeAlignatorDPFull( ALIGNMENT_LOCAL, 0, 0 ));
	HMultipleAlignator ma(makeMultipleAlignatorSimple( a ));
	std::vector<std::string>result1;
	result1.push_back("--------EEEEAAAADDDD--------MMMMEEEE--------");
	result1.push_back("----FFFF----AAAA----CCCC----MMMM----FFFF----");
	result1.push_back("GGGG--------AAAA--------KKKKMMMM--------GGGG");
	std::vector<std::string>result2;
	result2.push_back("----EEAAAADDDD--------MMMMEE----");
	result2.push_back("--FF--AAAA----CCCC----MMMM--FF--");
	result2.push_back("GG----AAAA--------KKKKMMMM----GG");

 	testMultipleAlignator( ma, result1, result2, false );
}

BOOST_AUTO_TEST_CASE( multiple_alignment_pileup )
{
	setDefaultSubstitutionMatrix( makeSubstitutionMatrix( getDefaultEncoder()->getAlphabetSize(), 1, -1000) );
	HAlignator a(makeAlignatorDPFull( ALIGNMENT_LOCAL, 0, 0 ));
	HMultipleAlignator ma(makeMultipleAlignatorPileup( a ));
	std::vector<std::string>result1;
	result1.push_back("EEEEAAAADDDDMMMMEEEE");
	result1.push_back("----AAAA----MMMM----");
	result1.push_back("----AAAA----MMMM----");
	std::vector<std::string>result2;
	result2.push_back("EEAAAADDDDMMMMEE");
	result2.push_back("--AAAA----MMMM--");
	result2.push_back("--AAAA----MMMM--");
	testMultipleAlignator( ma, result1, result2, true );
}
