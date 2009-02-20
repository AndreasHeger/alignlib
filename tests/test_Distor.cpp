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
#include "DistanceMatrix.h"
#include "HelpersDistanceMatrix.h"
#include "Distor.h"
#include "HelpersDistor.h"

using namespace std;
using namespace alignlib;

int main ()
{
  // create a multiple alignment
  HMultipleAlignment mali(makeMultipleAlignment());
  mali->add(alignlib::makeAlignatum("-AADDAACCAAA-"));
  mali->add(alignlib::makeAlignatum("AAKKAA-CCAAAA"));
  mali->add(alignlib::makeAlignatum("-A-AAA-CCA-A-"));
  mali->add(alignlib::makeAlignatum("AAAGAAA--AAAA"));

  HDistanceMatrix matrix(makeDistanceMatrixSymmetric(4, 0));

  // cout << *matrix << endl;

  HDistor d1(makeDistorKimura());
  d1->calculateMatrix( matrix, mali );
  // cout << *matrix << endl;

  HDistor d2(makeDistorClustal());
  d2->calculateMatrix( matrix, mali );
  // cout << *matrix << endl;

  return (EXIT_SUCCESS);

}
