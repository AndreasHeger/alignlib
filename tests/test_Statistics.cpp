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
#include <fstream>
#include <vector>

#include <time.h> 

#include "alignlib.h"

#include "Alignandum.h"
#include "HelpersSequence.h"

#include "Alignator.h"
#include "HelpersAlignator.h"

#include "Alignment.h"
#include "HelpersAlignment.h"

#include "Statistics.h"

#ifdef WITH_DMALLOC
#include <dmalloc.h>
#endif

using namespace std;

#ifdef WITH_DMALLOC
#include <dmalloc.h>
#endif

using namespace alignlib;

int main () {

  Alignator * a = makeFullDP( -10.0, -2.0);
  Alignandum * s1 = makeSequence( "AAACCCAAAAACCCAAAAAAA");
  Alignandum * s2 = makeSequence( "AAACCCAAAAACCCAAAAAAA");
    
  NormalDistributionParameters * result = makeNormalDistributionParameters();

  calculateZScoreParameters( result, s1, s2, a, 100);

  cout << "mean=" << result->getMean() << " std=" << result->getStandardDeviation() << endl;
  
  delete a;
  delete s1;
  delete s2;
  delete result;

}
