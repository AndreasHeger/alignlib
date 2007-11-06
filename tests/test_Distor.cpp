//--------------------------------------------------------------------------------
// Project LibAlign
//
// Copyright (C) 2000 Andreas Heger All rights reserved
//
// Author: Andreas Heger <heger@ebi.ac.uk>
//
// $Id: test_Distor.cpp,v 1.2 2004/06/02 12:14:35 aheger Exp $
//--------------------------------------------------------------------------------

/** Test Matrix object
*/

#if HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>

#include <time.h> 

#include "MultipleAlignment.h"
#include "HelpersMultipleAlignment.h"
#include "HelpersAlignatum.h"

#include "alignlib.h"
#include "PhyloMatrix.h"
#include "HelpersPhyloMatrix.h"
#include "Distor.h"
#include "HelpersDistor.h"

using namespace std;
using namespace alignlib;

#ifdef WITH_DMALLOC
#include <dmalloc.h>
#endif

int main ()
{
  // create a multiple alignment
  alignlib::MultipleAlignment * mali = alignlib::makeMultipleAlignment();
  mali->add(alignlib::makeAlignatumFromString("-AADDAACCAAA-"));
  mali->add(alignlib::makeAlignatumFromString("AAKKAA-CCAAAA"));
  mali->add(alignlib::makeAlignatumFromString("-A-AAA-CCA-A-"));
  mali->add(alignlib::makeAlignatumFromString("AAAGAAA--AAAA"));     

  PhyloMatrix * matrix = makePhyloMatrixSymmetric(4, 0);

  cout << *matrix << endl;

  Distor * d1 = makeDistorKimura();
  d1->calculateMatrix( matrix, mali );
  cout << *matrix << endl;

  Distor * d2 = makeDistorClustal();
  d2->calculateMatrix( matrix, mali );
  cout << *matrix << endl;

  delete d1;
  delete d2;
  delete matrix;
  delete mali;

  return (EXIT_SUCCESS);

}
