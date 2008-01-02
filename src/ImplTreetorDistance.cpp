//--------------------------------------------------------------------------------
// Project alignlib
//
// Copyright (C) 2000 Andreas Heger All rights reserved
//
// Author: Andreas Heger <heger@ebi.ac.uk>
//
// $Id: ImplTreetorDistance.cpp,v 1.1.1.1 2002/07/08 21:20:17 heger Exp $
//--------------------------------------------------------------------------------

#include <iostream>
#include <iomanip>
#include <cassert>
#include "MultipleAlignment.h"

#include "ImplTreetorDistance.h"
#include "Tree.h"
#include "PhyloMatrix.h"
#include "HelpersPhyloMatrix.h"
#include "Distor.h"
#include "AlignlibDebug.h"

using namespace std;

namespace alignlib {

//----------< constructors and destructors >--------------------------------------
ImplTreetorDistance::ImplTreetorDistance ( const Distor * distor ) : 
  ImplTreetor(),
  mDistor( distor ) {
}

ImplTreetorDistance::~ImplTreetorDistance () 
  {
}

ImplTreetorDistance::ImplTreetorDistance (const ImplTreetorDistance & src ) : ImplTreetor( src ) 
{
}

//------------------------------------------------------------------------------------------------------------------------------
void ImplTreetorDistance::startUp( Tree * tree, const alignlib::MultipleAlignment * mali) const 
{
	debug_func_cerr( 5 );

  // create a distance matrix
  if (mali != NULL) {
      assert( mali->getWidth() > 0);
      mWorkMatrix = makePhyloMatrixSymmetric( mali->getWidth());
  } else {
      mWorkMatrix = makePhyloMatrixSymmetric(0);
  }

  // call Distor for multiple alignment to fill matrix with distances
  mDistor->calculateMatrix( mWorkMatrix, mali );

  // save a pointer to the tree
  mTree = tree;
}

//------------------------------------------------------------------------------------------------------------------------------
void ImplTreetorDistance::cleanUp() const 
{

  /* I have to perform this check, as otherwise I got an error
     under gcc on freeing an empty pointer */
  if (mWorkMatrix != NULL)
    delete mWorkMatrix;
}

//------------------------------------------------------------------------------------------------------------------------------
void ImplTreetorDistance::swapHelpers( PhyloMatrixSize cluster_1, PhyloMatrixSize cluster_2) const 
{
}

//------------------------------------------------------------------------------------------------------------------------------
Tree * ImplTreetorDistance::calculateTree( Tree * tree, const alignlib::MultipleAlignment * mali ) const 
{

	debug_func_cerr( 5 );
  startUp( tree, mali );

  PhyloMatrixSize swap_temp;	// for swapping indices

#define SWAP(x,y) { swap_temp = x; x = y; y = swap_temp; }

  //------------------------------------------------------------------------------------------
  // start up 
  PhyloMatrixSize width = mWorkMatrix->getWidth();
  PhyloMatrixSize i;
  
  /* in this algorithm I assume that the matrix I use is only a half-matrix using the lower diagonal */
  tree->setNumLeaves( width );					// allocate memory for tree

  mIndices = new Node[width];			        /* allocate array to keep track of indices */

  for (i = 0; i < width; i++) mIndices[i] = i;
  
  //----------------------------------------------------------------------------------------------------
  // Perform hierarchical clustering

  PhyloMatrixSize last_row = width - 1;
  
  /* shrink distance matrix, until it contains only a single cluster */
  while (last_row > 0) 
  {
    
	  debug_cerr( 6, "Last row " << last_row );
	  debug_cerr( 6, "Work matrix " << endl << *mWorkMatrix );

    // find minimum distance in matrix 
    /*
	-
	--
	---
	----
	-x---
	------
	-------
	--------
      */
    calculateMinimumDistance();
    PhyloMatrixSize min_row = mMinimumCoord.row;
    PhyloMatrixSize min_col = mMinimumCoord.col;

    debug_cerr( 5, "Joining nodes -> "  
			    << "minimum distance :" << mMinimumValue << " " 
			    << "node 1: " << mIndices[min_row] << " (" << min_row << ") " 
			    << "node 2: " << mIndices[min_col] << " (" << min_col << ") " );
    
    //------------------------------------------------------------------------------------------------------
    // move rows around, so that the last two joined cluster are in the two last rows

    /*
	-			-
	y-			--
	-y-			---
	-y--		->	----
	xOxx-			-----
	-y--x-			------
	-y--x--			xxxxxxx
	-y--x---		yyyyyyyy
    */

    PhyloMatrixSize second_row = last_row - 1;		// second to last row

    // exchange row with last row
    mWorkMatrix->swap( min_row, last_row);
    SWAP( mIndices[min_row], mIndices[last_row]);
    swapHelpers( min_row, last_row);
    
    // exchange col with second to last row
    mWorkMatrix->swap( min_col, second_row);
    SWAP( mIndices[min_col], mIndices[second_row]);
    swapHelpers( min_col, second_row);

    //------------------------------------------------------------------------------------------------------
    // join the two nodes in the tree giving the correct edge weights
    Node new_node = joinNodes( second_row, last_row );
    
    debug_cerr( 6, "-> new node " << new_node ); 

    //------------------------------------------------------------------------------------------------------
    // calculate distance to new cluster and put them in second to last row
    updateDistanceMatrix( second_row, last_row);

    /* delete last row, update mIndices, so that the last row now
       contains the number of the new cluster id */
    /*
	-	       	-
	--	       	--
	---	       	---
	----	     ->	----
	-----		-----
	------		------
	xxxxxxx		xxxxxxx
	yyyyyyyy
    */
    
    mWorkMatrix->shrink();
    
    mIndices[last_row - 1] = new_node;
    
    last_row = mWorkMatrix->getWidth() - 1;
    
  }    
  
  delete [] mIndices;
  cleanUp();	// cleans up mWorkMatrix

  return mTree;
}

} /* namespace alignlib */
