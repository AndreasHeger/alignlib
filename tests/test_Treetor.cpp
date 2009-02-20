//--------------------------------------------------------------------------------
// Project LibAlign
//
// Copyright (C) 2000 Andreas Heger All rights reserved
//
// Author: Andreas Heger <heger@ebi.ac.uk>
//
// $Id: test_Treetor.cpp,v 1.2 2004/06/02 12:14:35 aheger Exp $
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

#include "alignlib.h"
#include "alignlib_fwd.h"

#include "MultipleAlignment.h"
#include "HelpersMultipleAlignment.h"
#include "HelpersAlignatum.h"
#include "DistanceMatrix.h"
#include "HelpersDistanceMatrix.h"
#include "Distor.h"
#include "HelpersDistor.h"
#include "Tree.h"
#include "HelpersTree.h"
#include "Treetor.h"
#include "HelpersTreetor.h"
#include "Matrix.h"

using namespace std;
using namespace alignlib;

DistanceMatrixValue source_linkage[]=  { 0.0, 4.0, 7.0, 4.0, 7.0,
				       4.0, 0.0, 7.0, 2.0, 7.0,
				       7.0, 7.0, 0.0, 7.0, 2.0,
				       4.0, 2.0, 7.0, 0.0, 7.0,
				       7.0, 7.0, 2.0, 7.0, 0.0 };


// the tree from the Graeme's lecture
/*  TYPE_MATRIX source_nj[] =  {
      0.0, 2.0, 4.0, 4.2,
      2.0, 0.0, 4.0, 4.0,
      4.0, 4.0, 0.0, 2.0,
      4.2, 4.0, 2.0, 0.0 };
*/


/** example from Durbin et al. p 170, result should be:
    ((2:0.4,0:0.1):0.05,(1:0.1,3:0.4):0.05)
*/

DistanceMatrixValue source_nj[] =  { 0.0, 0.3, 0.5, 0.6,
				   0.3, 0.0, 0.6, 0.5,
				   0.5, 0.6, 0.0, 0.9,
				   0.6, 0.5, 0.9, 0.0 };

void testTreetor( HTreetor & treetor )
{

  HTree tree (makeTree() );

  // create a multiple alignment
  HMultipleAlignment mali = makeMultipleAlignment();
  mali->add(makeAlignatum("-AADDAACCAAA-"));
  mali->add(makeAlignatum("AAKKAA-CCAAAA"));
  mali->add(makeAlignatum("-A-AAA-CCA-A-"));
  mali->add(makeAlignatum("AAAGAAA--AAAA"));

  treetor->calculateTree( tree, mali );
  // cout << *tree << endl;
  // writeNewHampshire( cout, tree);

}

int main ()
{

  /* test different tree building algorithms */

	HTreetor treetor;
	HDistor distor;
	HDistanceMatrix matrix;

  //------------------------> Test 1<-----------------------------------------
  // cout << "Test 1: create a tree from a multiple alignment:" << endl;
  treetor = makeTreetorDistanceLinkage( getDefaultDistor() );
  testTreetor( treetor );

  //------------------------> Test 2<-----------------------------------------
  // cout << "Test 2: creating a tree from a distance matrix:" << endl;
  matrix = makeDistanceMatrixSymmetric(5);
  fillDistanceMatrix( matrix, source_linkage );

  distor = makeDistorDummy( matrix );
  treetor = makeTreetorDistanceLinkage( distor, UPGMA );
  testTreetor( treetor );

  //------------------------> Test 2<-----------------------------------------
  // cout << "Test 3: creating a tree from a distance matrix:" << endl;
  matrix = makeDistanceMatrixSymmetric(4);
  fillDistanceMatrix( matrix, source_nj );

  distor = makeDistorDummy( matrix );
  treetor = makeTreetorDistanceNJ( distor );
  testTreetor( treetor );

  exit(EXIT_SUCCESS);

}
