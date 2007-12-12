//--------------------------------------------------------------------------------
// Project LibAlign
//
// Copyright (C) 2000 Andreas Heger All rights reserved
//
// Author: Andreas Heger <heger@ebi.ac.uk>
//
// $Id: ImplPhyloMatrixSymmetric.cpp,v 1.1.1.1 2002/07/08 21:20:17 heger Exp $
//--------------------------------------------------------------------------------    


#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdlib.h>
#include "alignlib.h"
#include "ImplPhyloMatrixSymmetric.h"
#include "AlignlibDebug.h"
#include "AlignException.h"

using namespace std;

namespace alignlib {

PhyloMatrix * makePhyloMatrixSymmetric( PhyloMatrixSize size, PhyloMatrixValue default_value)
{
	return new ImplPhyloMatrixSymmetric( size, default_value );
}

//---------------------------------------------------------< constructors and destructors >--------------------------------------

ImplPhyloMatrixSymmetric::ImplPhyloMatrixSymmetric() : ImplPhyloMatrix()
{
}

ImplPhyloMatrixSymmetric::~ImplPhyloMatrixSymmetric ()
{
}

// not calling the base class constructor is dangerous!!
ImplPhyloMatrixSymmetric::ImplPhyloMatrixSymmetric (PhyloMatrixSize width, PhyloMatrixValue default_value ) : ImplPhyloMatrix()
{
	debug_func_cerr( 5 );  

	mWidth  = width;
	mSize   = mWidth * (mWidth - 1) / 2;

#ifdef DEBUG
	cout << "Allocating " << mSize << " bytes for amatrix of width " << mWidth << endl;
#endif

	mMatrix = new PhyloMatrixValue[ mSize ];

	PhyloMatrixSize i;
	for (i = 0; i < mSize; i++)
		mMatrix[i] = default_value;
}

ImplPhyloMatrixSymmetric::ImplPhyloMatrixSymmetric( const ImplPhyloMatrixSymmetric & src) : ImplPhyloMatrix( src )
{

	debug_func_cerr( 5 );

	//!! to be tightened up. Do not copy twice, otherwise memory leak!!!
	// since virtual functions do not work in constructors, I have to everything myself

	mWidth = src.getWidth();
	const PhyloMatrixValue * matrix = src.mMatrix;

	mSize = mWidth * (mWidth - 1) / 2;
	mMatrix = new PhyloMatrixValue[ mSize ];

	if (!mMatrix)
		throw AlignException("Out of memory in ImplPhyloMatrixSymmetric");

	PhyloMatrixSize i, j;
	PhyloMatrixSize index = 0;

	for (i = 1; i < mWidth; i++)			// iterate through rows
		for (j = 0; j < i; j++) 			// iterate through columns
			mMatrix[index++] = matrix[src.getIndex(i,j)];
}

//-------------------------------------------------------------
void ImplPhyloMatrixSymmetric::shrink() 
{
	debug_func_cerr( 5 );
	/* simply decrease the width of the matrix */
	mWidth --;
	mSize -= mWidth;

	debug_cerr( 5, "New matrix size: " << mWidth << "(" << mSize << ")" );

}

//-------------------------------------------------------------
void ImplPhyloMatrixSymmetric::swap( PhyloMatrixSize row_1, PhyloMatrixSize row_2 )
{
	debug_func_cerr( 5 );

	/*
	-			-
	y-			x-
	-y-			-x-
	-y--		->	-x--
	xOxx-			yyyy-
	-y--x-			-x--y-
	-y--x--			-x--y--
	-y--x---		-x--y---
   i.e. two columns have to skipped:
   1. col = row_1
   2. col = row_2
	 */

	PhyloMatrixValue t;

#define SWAP( x, y) { t = mMatrix[x], mMatrix[x] = mMatrix[y], mMatrix[y] = t;}

	PhyloMatrixSize i;

	if (row_2 < row_1) 
	{
		i = row_1; row_1 = row_2; row_2 = i;
	}

	for (i = 0; i < row_1; i++)
		SWAP( getIndex(row_1,i), getIndex( row_2,i));

	for (i = row_1 + 1; i < row_2; i++) 
		SWAP( getIndex(row_1,i), getIndex( row_2,i));

	for (i = row_2 + 1; i < mWidth; i++)
		SWAP( getIndex(row_1,i), getIndex( row_2,i));

#undef SWAP

}

//-------------------------------------------------------------
PhyloMatrixSize ImplPhyloMatrixSymmetric::getIndex( PhyloMatrixSize row, PhyloMatrixSize col) const
{
	debug_func_cerr( 5 );

	int x;

	if (row == col)
		return NO_INDEX;           // the diagonal is not part of the matrix

	if (row > col)
		x = row * ( row - 1) / 2 + col;
	else
		x = col * ( col - 1) / 2 + row;

	return x;
}

//-------------------------------------------------------------
PhyloMatrixSize ImplPhyloMatrixSymmetric::getRow( PhyloMatrixSize index ) const
{
	debug_func_cerr( 5 );

	PhyloMatrixSize row = 1;
	PhyloMatrixSize row_index = 0;

	while (index >= row_index)
		row_index += row++;

	return (row - 1);
}


//-------------------------------------------------------------
PhyloMatrixSize ImplPhyloMatrixSymmetric::getColumn( PhyloMatrixSize index ) const {
	// computationally inefficient, maybe use cache
	return (index - getIndex(getRow(index), 0));	
}

//-------------------------------------------------------------
void ImplPhyloMatrixSymmetric::calculateSize()
{

	debug_func_cerr( 5 );

	mSize = mWidth * (mWidth - 1) / 2;
}

//--------------------------------------------------------------------------------------------------------------------------------


} // namespace alignlib
