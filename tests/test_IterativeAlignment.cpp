/*
  alignlib - a library for aligning protein sequences

  $Id: test_AlignatorDots.cpp,v 1.3 2004/06/02 12:11:38 aheger Exp $

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
#include <vector>

#include <time.h> 

#include "alignlib.h"


#include "HelpersSequence.h"

#include "HelpersProfile.h"

#include "Alignator.h"
#include "HelpersAlignator.h"

#include "Alignata.h"
#include "HelpersAlignata.h"

#ifdef WITH_DMALLOC
#include <dmalloc.h>
#endif

using namespace std;

#ifdef WITH_DMALLOC
#include <dmalloc.h>
#endif

using namespace alignlib;

void test() 
{
  {
	  cout << "---------------------Testing Alignator----------------------------------" << endl;
	  // build alignment

	  Alignandum * row = makeSequence( "AAACCCCAAACCCCCCCAAACCCCCCCAAACCCCCCCAAA" );
      Alignandum * col = makeSequence( "AAAKKKKAAAKKKKKKKAAAKKKKKKKAAAKKKKKKKAAA" );

      Alignata * result = makeAlignataVector();

      Alignator * alignator = makeAlignatorDPFull( ALIGNMENT_LOCAL, -10.0, -2.0 );
      alignator->align( row, col, result);

      performIterativeAlignment( result, row, col,
    		  alignator, 1.0);

  }
}


int main () {

  test();

  return EXIT_SUCCESS;
}
