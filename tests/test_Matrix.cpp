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
#include <vector>

#include "Matrix.h"

#define BOOST_TEST_MODULE
#include <boost/test/included/unit_test.hpp>
using boost::unit_test::test_suite;

using namespace std;
using namespace alignlib;

BOOST_AUTO_TEST_CASE( test_map1 )
{
	Matrix<int>matrix( 3, 3, 1);
	std::cout << matrix << std::endl;
	std::vector<unsigned int> map_new2old( 2 );
	map_new2old[0] = 0;
	map_new2old[1] = 2;
	matrix.permuteRows( map_new2old );
	std::cout << matrix << std::endl;
	matrix.permuteCols( map_new2old );
	std::cout << matrix << std::endl;	
}

BOOST_AUTO_TEST_CASE( test_map2 )
{
	Matrix<int>matrix( 2, 2, 1);
	std::cout << matrix << std::endl;
	std::vector<unsigned int> map_new2old( 2 );
	map_new2old[0] = 0;
	map_new2old[1] = 1;
	matrix.permuteRows( map_new2old );
	std::cout << matrix << std::endl;
	matrix.permuteCols( map_new2old );
	std::cout << matrix << std::endl;	
}










