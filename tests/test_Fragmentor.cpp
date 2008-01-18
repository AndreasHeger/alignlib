/*
  alignlib - a library for aligning protein sequences

  $Id: test_Fragmentor.cpp,v 1.3 2004/06/02 12:11:38 aheger Exp $

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
#include "alignlib_fwd.h"

#include "HelpersSequence.h"
#include "HelpersProfile.h"
#include "Alignandum.h"
#include "Alignator.h"
#include "HelpersAlignator.h"
#include "Fragmentor.h"
#include "HelpersFragmentor.h"
#include "Alignment.h"
#include "HelpersAlignment.h"
#include "HelpersEncoder.h"

using namespace std;
using namespace alignlib;

int main () 
{

  HAlignment dots(makeAlignmentMatrixRow());

  HAlignandum s2(makeSequence("AAAA", getDefaultEncoder() ) );
  HAlignandum s1(makeSequence("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA", 
		  getDefaultEncoder() ) );

  //  std::ifstream fin("test.dots");
  // readAlignmentPairs( dots, fin );

  fillAlignmentIdentity( dots, 1, 4, 0);
  fillAlignmentIdentity( dots, 10, 14, 0);
  fillAlignmentIdentity( dots, 20, 24, 0);
  rescoreAlignment( dots, 1.0);

  cout << "Dots read in:\n" << *dots << endl;

  HFragmentor f(makeFragmentorIterative( dots, 2, -1, -1));
  
  HAlignment t_ali(makeAlignmentSet());
    
  HFragmentVector fragments(f->fragment( t_ali, s1, s2 ));

  for (unsigned int i = 0; i < fragments->size(); i++) {
    cout << "fragment " << i << endl << *((*fragments)[i]);
  }
  
}
