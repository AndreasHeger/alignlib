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

BOOST_AUTO_TEST_CASE( multiple_alignment_simple )
{
	setDefaultSubstitutionMatrix( makeSubstitutionMatrix( getDefaultEncoder()->getAlphabetSize(), 1, -1000) );
	HAlignator a(makeAlignatorDPFull( ALIGNMENT_LOCAL, 0, 0 ));
	HStringVector sequences(new StringVector());
	sequences->push_back( "EEEEAAAADDDDMMMMEEEE");
	sequences->push_back( "FFFFAAAACCCCMMMMFFFF");
	sequences->push_back( "GGGGAAAAKKKKMMMMGGGG");
	HMultAlignment result(makeMultAlignment());

	HMultipleAlignator ma(makeMultipleAlignatorSimple( a ));
	ma->align( result, sequences );

	MultAlignmentFormatPlain f( result, sequences);
	BOOST_CHECK_EQUAL( f.mData[0]->getString(), "--------EEEEAAAADDDD--------MMMMEEEE--------");
	BOOST_CHECK_EQUAL( f.mData[1]->getString(), "----FFFF----AAAA----CCCC----MMMM----FFFF----");
	BOOST_CHECK_EQUAL( f.mData[2]->getString(), "GGGG--------AAAA--------KKKKMMMM--------GGGG");
}

BOOST_AUTO_TEST_CASE( multiple_alignment_pileup )
{
	setDefaultSubstitutionMatrix( makeSubstitutionMatrix( getDefaultEncoder()->getAlphabetSize(), 1, -1000) );
	HAlignator a(makeAlignatorDPFull( ALIGNMENT_LOCAL, 0, 0 ));
	HStringVector sequences(new StringVector());
	sequences->push_back( "EEEEAAAADDDDMMMMEEEE");
	sequences->push_back( "FFFFAAAACCCCMMMMFFFF");
	sequences->push_back( "GGGGAAAAKKKKMMMMGGGG");
	HMultAlignment result(makeMultAlignment());

	HMultipleAlignator ma(makeMultipleAlignatorPileup( a ));
	ma->align( result, sequences );

	MultAlignmentFormatPlain f( result, sequences);
	BOOST_CHECK_EQUAL( f.mData[0]->getString(), "EEEEAAAADDDDMMMMEEEE");
	BOOST_CHECK_EQUAL( f.mData[1]->getString(), "----AAAA----MMMM----");
	BOOST_CHECK_EQUAL( f.mData[2]->getString(), "----AAAA----MMMM----");
}

