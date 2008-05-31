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
#include "AlignmentFormat.h"
#include "HelpersAlignment.h"
#include "AlignlibDebug.h"


#define BOOST_TEST_MODULE
#include <boost/test/included/unit_test.hpp>
#include <boost/test/unit_test_monitor.hpp>
using boost::unit_test::test_suite;
/*
#include <boost/test/unit_test.hpp>
using namespace boost::unit_test;
*/
using namespace std;
using namespace alignlib;


void throw_AlignlibException()
{
    throw AlignlibException( "AlignlibException test" );
}

void AlignlibException_translator( AlignlibException )
{
    BOOST_TEST_MESSAGE( "Caught Alignlib Exception" );
}

BOOST_AUTO_TEST_CASE( test_AlignlibException )
{
    boost::unit_test::unit_test_monitor.register_exception_translator<AlignlibException>( &AlignlibException_translator );
    BOOST_TEST_CASE( &throw_AlignlibException );
}
