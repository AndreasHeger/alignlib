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

#include "MultipleAlignment.h"
#include "HelpersMultipleAlignment.h"
#include "HelpersAlignatum.h"
#include "PhyloMatrix.h"
#include "HelpersPhyloMatrix.h"
#include "Distor.h"
#include "HelpersDistor.h"
#include "Tree.h"
#include "HelpersTree.h"
#include "Treetor.h"
#include "HelpersTreetor.h"

#ifdef WITH_DMALLOC
#include <dmalloc.h>
#endif

using namespace std;
using namespace alignlib;

PhyloMatrixValue source_linkage[]=  { 0.0, 4.0, 7.0, 4.0, 7.0,
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

PhyloMatrixValue source_nj[] =  { 0.0, 0.3, 0.5, 0.6,
				   0.3, 0.0, 0.6, 0.5,
				   0.5, 0.6, 0.0, 0.9,
				   0.6, 0.5, 0.9, 0.0 };

void testTreetor( Treetor * treetor ) {

  Tree * tree = makeTree();
  
  // create a multiple alignment
  alignlib::MultipleAlignment * mali = alignlib::makeMultipleAlignment();
  mali->add(alignlib::makeAlignatumFromString("-AADDAACCAAA-"));
  mali->add(alignlib::makeAlignatumFromString("AAKKAA-CCAAAA"));
  mali->add(alignlib::makeAlignatumFromString("-A-AAA-CCA-A-"));
  mali->add(alignlib::makeAlignatumFromString("AAAGAAA--AAAA"));     

  treetor->calculateTree( tree, mali );
  cout << *tree << endl;
  // writeNewHampshire( cout, tree);
  delete mali;

}

int main () {

  /* test different tree building algorithms */
  
  Treetor * treetor;
  Distor * distor;
  PhyloMatrix * matrix;

  //------------------------> Test 1<-----------------------------------------
  cout << "Test 1: create a tree from a multiple alignment:" << endl;
  treetor = makeTreetorDistanceLinkage();
  testTreetor( treetor );
  delete treetor;
  
  //------------------------> Test 2<-----------------------------------------
  cout << "Test 2: creating a tree from a distance matrix:" << endl;
  matrix = makePhyloMatrixSymmetric(5);
  fillPhyloMatrix( matrix, source_linkage );
  
  distor = makeDistorDummy( matrix );
  treetor = makeTreetorDistanceLinkage( UPGMA, distor );
  testTreetor( treetor );
  delete distor;
  delete matrix;
  delete treetor;
  
  //------------------------> Test 2<-----------------------------------------
  cout << "Test 3: creating a tree from a distance matrix:" << endl;
  matrix = makePhyloMatrixSymmetric(4);
  fillPhyloMatrix( matrix, source_nj );
  
  distor = makeDistorDummy( matrix );
  treetor = makeTreetorDistanceNJ( distor );
  testTreetor( treetor );
  delete distor;
  delete matrix;
  delete treetor;

  exit(EXIT_SUCCESS);

}
