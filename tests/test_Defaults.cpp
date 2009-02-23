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

#define DEFAULT_TEST( test, handle, get ) \
BOOST_AUTO_TEST_CASE( test ) \
{ \
	const handle & l1 = getDefaultToolkit()->get(); \
	const handle & l2 = l1->getToolkit()->get(); \
	BOOST_CHECK_EQUAL( l1, l2); \
}

DEFAULT_TEST( test_Weightor, HWeightor, getWeightor );
DEFAULT_TEST( test_Encoder, HEncoder, getEncoder );
DEFAULT_TEST( test_Alignator, HAlignator, getAlignator );
DEFAULT_TEST( test_Alignment, HAlignment, getAlignment );
DEFAULT_TEST( test_MultAlignment, HMultAlignment, getMultAlignment );
DEFAULT_TEST( test_MultipleAlignator, HMultipleAlignator, getMultipleAlignator );
DEFAULT_TEST( test_LogOddor, HLogOddor, getLogOddor );
DEFAULT_TEST( test_Regularizor, HRegularizor, getRegularizor );
DEFAULT_TEST( test_Distor, HDistor, getDistor );
DEFAULT_TEST( test_Treetor, HTreetor, getTreetor );
DEFAULT_TEST( test_Fragmentor, HFragmentor, getFragmentor );
DEFAULT_TEST( test_Scorer, HScorer, getScorer);
DEFAULT_TEST( test_Iterator2D, HIterator2D, getIterator2D );
DEFAULT_TEST( test_SubstitutionMatrix, HSubstitutionMatrix, getSubstitutionMatrix );



