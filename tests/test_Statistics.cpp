/*
  alignlib - a library for aligning protein sequences

  $Id: test_Statistics.cpp,v 1.3 2004/06/02 12:11:38 aheger Exp $

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

#include "alignlib.h"

using namespace std;
using namespace alignlib;

#define BOOST_TEST_MODULE
#include <boost/test/included/unit_test.hpp>
using boost::unit_test::test_suite;

BOOST_AUTO_TEST_CASE( test_LogOddorUniform )
{
	HAlignator a = makeAlignatorDPFull( 
			ALIGNMENT_GLOBAL,
			-10.0, -2.0);
	HAlignandum s1 = makeSequence( "AAACCCAAAAACCCAAAAAAA");
	HAlignandum s2 = makeSequence( "AAACCCAAAAACCCAAAAAAA");
	    
	NormalDistributionParameters * result = makeNormalDistributionParameters();

	calculateZScoreParameters( result, s1, s2, a, 100);
	cout << "mean=" << result->getMean() << " std=" << result->getStandardDeviation() << endl;
}

