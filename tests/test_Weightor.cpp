/*
  alignlib - a library for aligning protein sequences

  $Id: test_MultipleAlignment.cpp,v 1.6 2004/06/02 12:11:38 aheger Exp $

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

/** Test the MultipleAlignment - object
 */

#if HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>
#include <fstream>
#include <sstream>
#include <cassert>

#include "alignlib.h"
#include "alignlib_fwd.h"

#include "Alignatum.h"
#include "Alignandum.h"
#include "HelpersAlignatum.h"
#include "MultipleAlignment.h"
#include "HelpersMultipleAlignment.h"
#include "AlignlibDebug.h"
#include "Weightor.h"
#include "HelpersEncoder.h"

#define BOOST_TEST_MODULE
#include <boost/test/included/unit_test.hpp>
using boost::unit_test::test_suite;

using namespace std;
using namespace alignlib;

std::string ref_protein20 = "ACDEFGHIKLMNPQRSTVWY";
std::string ref_protein20x3 = ref_protein20 + ref_protein20 + ref_protein20;

void test_GenericWeightor( HWeightor & weightor )
{
	HAlignandum a(makeProfile( ref_protein20x3, 3,
			getDefaultEncoder(),
			weightor,
			getDefaultRegularizor(),
			getDefaultLogOddor()));

	a->prepare();
}

BOOST_AUTO_TEST_CASE( test_Weightor )
{
	HWeightor l(makeWeightor());
	test_GenericWeightor( l );
}

BOOST_AUTO_TEST_CASE( test_WeightorHenikoff )
{
	HWeightor l(makeWeightorHenikoff());
	test_GenericWeightor( l );
}






