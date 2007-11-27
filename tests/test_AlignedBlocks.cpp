/*
  alignlib - a library for aligning protein sequences

  $Id: test_HelpersAlignata.cpp,v 1.5 2004/10/14 23:34:09 aheger Exp $

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

#include "alignlib.h"

#include "Alignata.h"
#include "HelpersAlignata.h"
#include "AlignataIterator.h"
#include "AlignedBlocks.h"

using namespace std;

using namespace alignlib;

bool isIdentical( const Alignata * a, const Alignata * b, bool inverse = false) 
{

  AlignataConstIterator it1(a->begin());
  AlignataConstIterator it1_end(a->end());

  AlignataConstIterator it2(b->begin());
  AlignataConstIterator it2_end(b->end());

  bool is_identical = true;

  for (; it1 != it1_end; ++it1, ++it2) {
    if (!inverse) {
      if (it1->mRow != it2->mRow && it1->mCol != it2->mCol) 
        is_identical = false;
    } else {
      if (it1->mRow != it2->mCol && it1->mCol != it2->mRow) 
        is_identical = false;
    }
  }

  return is_identical;
}

int main () 
{

	cout << "---------------------Testing fill-------------------------------" << endl;
	Alignata * a = makeAlignataVector();
	fillAlignataIdentity( a, 5, 10, 0);
	fillAlignataIdentity( a, 10, 15, 5);
	fillAlignataIdentity( a, 25, 30, -5);
	
	AlignedBlocks blocks_out( a );
	
	stringstream stream( stringstream::in | stringstream::out);
	stream << blocks_out;
	
	AlignedBlocks blocks_in;
	stream >> blocks_in;
	
	Alignata * b = makeAlignataVector();
	blocks_in.copy( b );

	assert( isIdentical( a, b) );
	
	delete a;
	delete b;
	
	exit(EXIT_SUCCESS);
}








