//--------------------------------------------------------------------------------
// Project LibAlign
//
// Copyright (C) 2000 Andreas Heger All rights reserved
//
// Author: Andreas Heger <heger@ebi.ac.uk>
//
// $Id: ImplPhyloMatrix.cpp,v 1.1.1.1 2002/07/08 21:20:17 heger Exp $
//--------------------------------------------------------------------------------    


#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdlib.h>
#include <cstring> // for memcpy
#include <cassert>
#include "alignlib_fwd.h"
#include "alignlib_interfaces.h"
#include "ImplPhyloMatrix.h"
#include "AlignlibDebug.h"
#include "AlignException.h"

using namespace std;

namespace alignlib {

//define swap_temp for this macro to work
#define SWAP(x,y) { swap_temp = x; x = y; y = swap_temp; }

//---------------------------------------------------------< constructors and destructors >--------------------------------------

ImplPhyloMatrix::ImplPhyloMatrix() : mWidth(0), mSize(0), mMatrix( NULL)
{
	debug_func_cerr( 5 );
}

ImplPhyloMatrix::ImplPhyloMatrix( PhyloMatrixSize width, PhyloMatrixValue default_value)
{
	debug_func_cerr( 5 );

	allocateMemory();

	PhyloMatrixSize i;
	for (i = 0; i < mSize; i++) 
		mMatrix[i] = default_value;
}

ImplPhyloMatrix::ImplPhyloMatrix( PhyloMatrixSize width, PhyloMatrixValue * source)
{
	debug_func_cerr( 5 );

	mWidth  = width;
	mMatrix = source;

	calculateSize();
}


ImplPhyloMatrix::~ImplPhyloMatrix ()
{
	debug_func_cerr( 5 );

	freeMemory();
}


ImplPhyloMatrix::ImplPhyloMatrix (const ImplPhyloMatrix & src ) : 
	mWidth( src.mWidth ), mSize( src.mSize )
	{

	debug_func_cerr( 5 );

	// have to allocate memory myself, because virtual functions do not work inside constructors
	// so calling CalculateSize() is not an option

	mMatrix = new PhyloMatrixValue[ mSize ];

	if (!mMatrix)
		throw AlignException("Out of memory in ImplPhyloMatrix::AllocateMemory");

	memcpy( mMatrix, src.mMatrix, mSize * sizeof( PhyloMatrixValue) );

	}


//--------------------------------------------------------------------------------------------
PhyloMatrixValue ImplPhyloMatrix::getMinimum( Coordinate & x) const
{
	debug_func_cerr( 5 );

	PhyloMatrixValue min = std::numeric_limits<PhyloMatrixValue>::max();

	PhyloMatrixSize i;
	PhyloMatrixSize best_index = 0;

	for (i = 0; i < mSize; i++) 
	{
		if (mMatrix[i] < min) 
		{
			min = mMatrix[i];
			best_index = i;
		}
	}

	x.row = getRow( best_index );
	x.col = getColumn( best_index );

	return min;
}

//--------------------------------------------------------------------------------------------
PhyloMatrixValue ImplPhyloMatrix::getMinimum() const
{
	debug_func_cerr( 5 );

	Coordinate x;

	return getMinimum(x);
}

//-------------------------------------------------------------
void ImplPhyloMatrix::setWidth(PhyloMatrixSize width)
{
	debug_func_cerr( 5 );

	mWidth = width;
	allocateMemory();
}

//--------------------------------------------------------------------------------------------
PhyloMatrixValue ImplPhyloMatrix::getMaximum( Coordinate & x) const
{ 
	debug_func_cerr( 5 );

	PhyloMatrixValue max = -999999;

	PhyloMatrixSize i;
	PhyloMatrixSize best_index = 0;

	for (i = 0; i < mSize; i++) {
		if (mMatrix[i] > max)
		{
			max = mMatrix[i];
			best_index = i;
		}
	}

	x.row = getRow( best_index );
	x.col = getColumn( best_index );
	return max;
}

PhyloMatrixValue ImplPhyloMatrix::getMaximum() const 
{
	Coordinate x;

	return getMaximum(x);
}


//--------------------------------------------------------------------------------------------
void ImplPhyloMatrix::allocateMemory()
{
	debug_func_cerr( 5 );

	// clear old memory, be careful!, use copy constructor to copy matrices
	PhyloMatrixSize saved_width = mWidth;
	freeMemory();
	mWidth = saved_width;

	debug_cerr( 5, "Allocating " << mSize << " bytes for amatrix of width " << mWidth );

	calculateSize();
	mMatrix = new PhyloMatrixValue[ mSize ];

	if (!mMatrix)
		throw AlignException("Out of memory in ImplPhyloMatrix::allocateMemory");

}

//--------------------------------------------------------------------------------------------
void ImplPhyloMatrix::freeMemory()
{
	debug_func_cerr( 5 );

	if (mMatrix != NULL) 
		delete [] mMatrix;

	mMatrix = NULL;
	mWidth  = 0;
	mSize   = 0;

}

//-------------------------------------------------------------
PhyloMatrixSize ImplPhyloMatrix::getWidth() const 
{
	return mWidth;
}

//-------------------------------------------------------------
PhyloMatrixSize ImplPhyloMatrix::getSize() const 
{
	return mSize;
}

//-------------------------------------------------------------
PhyloMatrixValue ImplPhyloMatrix::operator()( PhyloMatrixSize row, PhyloMatrixSize col) const 
{

	debug_func_cerr( 5 );

	debug_cerr( 5, "Looking up " << row << " " << col << ": index is " << getIndex( row, col ) ); 

	return mMatrix[getIndex( row, col)];
}    

//-------------------------------------------------------------
PhyloMatrixValue ImplPhyloMatrix::getElement( PhyloMatrixSize row, PhyloMatrixSize col) const 
{
	return mMatrix[getIndex( row, col)];
}    

//-------------------------------------------------------------
PhyloMatrixValue & ImplPhyloMatrix::operator()( PhyloMatrixSize row, PhyloMatrixSize col) 
{
	return mMatrix[getIndex(row,col)];
}    
//-------------------------------------------------------------
void ImplPhyloMatrix::setElement( PhyloMatrixSize row, PhyloMatrixSize col, PhyloMatrixValue value) 
{
	mMatrix[getIndex(row,col)] = value;
}    

//-------------------------------------------------------------
PhyloMatrixSize ImplPhyloMatrix::getIndex( PhyloMatrixSize row, PhyloMatrixSize col) const 
{
	return (row * mWidth + col);
}

//-------------------------------------------------------------
PhyloMatrixSize ImplPhyloMatrix::getRow( PhyloMatrixSize index ) const 
{
	return ( (PhyloMatrixSize)(index / mWidth)  ); 
}

//-------------------------------------------------------------
PhyloMatrixSize ImplPhyloMatrix::getColumn( PhyloMatrixSize index ) const 
{
	return ( index % mWidth );
}

//-------------------------------------------------------------
void ImplPhyloMatrix::calculateSize() 
{

	debug_func_cerr(5);

	mSize = mWidth * mWidth;
}

//---------------------------------------------------------< Input/Output routines >---------------------------------------------
void ImplPhyloMatrix::write( std::ostream & output ) const 
{
	// write so that output can be used in phylip

	cout << " " << mWidth << " " << mSize << endl;

	PhyloMatrixSize i, j;

	for (i = 0; i < mWidth; i++) {
		cout << i << "\t";
		for (j = 0; j < mWidth; j++) 
			cout << setw(10) << setprecision(4) << operator()( i, j) << " ";
		cout << endl;
	}
}

//---------------------------------------------------------< Input/Output routines >---------------------------------------------
void ImplPhyloMatrix::read( std::istream & input ) const 
{
}

//--------------------------------------------------------------------------------------------------------------------------------


} // namespace alignlib
