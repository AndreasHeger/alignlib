//--------------------------------------------------------------------------------
// Project alignlib
//
// Copyright (C) 2000 Andreas Heger All rights reserved
//
// Author: Andreas Heger <heger@ebi.ac.uk>
//
// $Id: ImplTreetorDistanceNJ.cpp,v 1.1.1.1 2002/07/08 21:20:17 heger Exp $
//--------------------------------------------------------------------------------

#include <iostream>
#include <iomanip>
#include <cassert>
#include "ImplTreetorDistanceNJ.h"
#include "Tree.h"
#include "PhyloMatrix.h"
#include "Distor.h"
#include "HelpersDistor.h"
#include "AlignException.h"
#include "AlignlibDebug.h"

using namespace std;

namespace alignlib {

//------------------------------------------------------< factory functions >------------------------------------------------
Treetor * makeTreetorDistanceNJ( const Distor * distor ) 
{
	if (distor == NULL) 
		return new ImplTreetorDistanceNJ( getDefaultDistor());
	else 
		return new ImplTreetorDistanceNJ( distor );
}

//---------------------------------------------------------< constructors and destructors >--------------------------------------
ImplTreetorDistanceNJ::ImplTreetorDistanceNJ ( const Distor * distor) : 
	ImplTreetorDistance( distor), mR (NULL) 
{
}

ImplTreetorDistanceNJ::~ImplTreetorDistanceNJ () 
{
}

ImplTreetorDistanceNJ::ImplTreetorDistanceNJ (const ImplTreetorDistanceNJ & src ) : 
	ImplTreetorDistance(src), mR( NULL) 
{
	// should not be able to copy while object is performing action.
	assert( src.mR == NULL );
}

//------------------------------------------------------------------------------------------------------------------------------
void ImplTreetorDistanceNJ::startUp( Tree * tree, const alignlib::MultipleAlignment * mali) const 
{

	cleanUp();
	ImplTreetorDistance::startUp( tree, mali);

	PhyloMatrixSize width = mWorkMatrix->getWidth();
	PhyloMatrixSize width_2 = width - 2;

	mR = new PhyloMatrixValue[ width ];

	PhyloMatrixSize i,k;

	for (i = 0; i < width; i++) {
		mR[i] = 0;
		for (k = 0; k < width; k++) 
			mR[i] += (*mWorkMatrix)( i, k );
		mR[i] /= width_2;
	}
}

//------------------------------------------------------------------------------------------------------------------------------
void ImplTreetorDistanceNJ::cleanUp() const 
{
	if (mR != NULL) 
		delete [] mR;

	mR = NULL;
}


//------------------------------------------------------------------------------------------------------------------------------
Node ImplTreetorDistanceNJ::joinNodes( PhyloMatrixSize cluster_1, PhyloMatrixSize cluster_2 ) const 
{
	debug_func_cerr( 5 );	

	PhyloMatrixValue d_ij = (*mWorkMatrix)( cluster_1, cluster_2 );
	PhyloMatrixValue d_ik = (d_ij + mR[cluster_1] - mR[cluster_2]) / 2;
	PhyloMatrixValue d_jk = d_ij - d_ik;

	debug_cerr( 5, "Joining nodes " <<
			mIndices[cluster_1] << " (" << d_ik << ") with " << 
			mIndices[cluster_2] << " (" << d_jk << ");" );

	return mTree->joinNodes( mIndices[ cluster_1 ], 
			mIndices[ cluster_2 ],
			d_ik, 
			d_jk );
};

//------------------------------------------------------------------------------------------------------------------------------
void ImplTreetorDistanceNJ::swapHelpers( PhyloMatrixSize cluster_1, PhyloMatrixSize cluster_2) const 
{
	PhyloMatrixValue t;

	t = mR[cluster_1];
	mR[cluster_1] = mR[cluster_2];
	mR[cluster_2] = t;

}

//------------------------------------------------------------------------------------------------------------------------------
void ImplTreetorDistanceNJ::calculateMinimumDistance() const 
{

	PhyloMatrixValue min = 999999;

	PhyloMatrixSize row, col;
	PhyloMatrixSize best_row = 0;
	PhyloMatrixSize best_col = 0;
	PhyloMatrixSize width = mWorkMatrix->getWidth();
	PhyloMatrixValue d; 

	// iterate through rows and columns, not through indices, since ri and rj have
	// to be calculated, and looking up the row/col for a given index is slow.

	for (row = 0; row < width - 1; row++) {
		for (col = row + 1; col < width; col++) {
			if ( (d = (*mWorkMatrix)( row, col ) - (mR[row] + mR[col])) < min ) {
				min = d;
				best_row = row;
				best_col = col;
			}
		}
	}

	mMinimumCoord.row = best_row;
	mMinimumCoord.col = best_col;

	mMinimumValue = min;
}

//------------------------------------------------------------------------------------------------------------------------------
void ImplTreetorDistanceNJ::updateDistanceMatrix( PhyloMatrixSize cluster_1, PhyloMatrixSize cluster_2 ) const 
{

	debug_func_cerr( 5 );
	// calculate distance to new cluster and put them in cluster_1

	//------------------------------------------------------------------------------------------------------
	// see Durbin et al. for description of algorithm and symbols used here.
	// the update proceeds in three steps
	// 1. mR[s]: multiply by |L| - 2 and subract d_si and d_sj
	// 2. calculate new distance d_sk (k = new cluster, is stored in cluster_1, i.e. i)
	// 3. mR[s]: add d_sk , divide by |L| - 3. 

	PhyloMatrixSize s;

	PhyloMatrixValue new_r = 0;		// mR for new cluster

	PhyloMatrixSize last_row = mWorkMatrix->getWidth() - 1;
	PhyloMatrixSize num_leaves_2 = last_row - 1;
	PhyloMatrixSize num_leaves_3 = last_row - 2;

	if (num_leaves_3 == 0)
		num_leaves_3 = 1;			

	// iterate to lastrow. if s = cluster_1 d_ss, etc are 0, so nothing changes.
	PhyloMatrixValue d_ij = (*mWorkMatrix)( cluster_1, cluster_2);

	for (s = 0; s < last_row; s++) {

		// update mR
		mR[s] =  mR[s] * num_leaves_2 
		- (*mWorkMatrix)( cluster_1, s) 
		- (*mWorkMatrix)( cluster_2, s);

		// calculate new value for d_ms
		PhyloMatrixValue d_is = (*mWorkMatrix)( cluster_1, s );
		PhyloMatrixValue d_js = (*mWorkMatrix)( cluster_2, s );

		PhyloMatrixValue new_dist = (d_is + d_js - d_ij) / 2.0;

		(*mWorkMatrix)( cluster_1, s) = new_dist; 
		new_r += new_dist;

		// update mR
		mR[s] = (mR[s] +  new_dist) / num_leaves_3;
	}    

	mR[cluster_1] = new_r / num_leaves_3;


}


} /* namespace alignlib */

