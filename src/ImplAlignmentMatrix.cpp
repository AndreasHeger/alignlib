/*
  alignlib - a library for aligning protein sequences

  $Id: ImplAlignmentMatrix.cpp,v 1.6 2004/09/24 19:03:27 aheger Exp $

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


#include <iostream> 
#include <algorithm>
#include <set>
#include <cassert>
#include "alignlib.h"
#include "AlignlibDebug.h"
#include "ImplAlignmentMatrix.h"
#include "AlignmentIterator.h"
#include "AlignException.h"

using namespace std;

namespace alignlib
{

/** Note: There are two arrays that have to be taken care of: mPairs and mIndex. If you
      need to copy one, you have to copy the other one as well. The same applies to deletion.

      clear() does not free the memory for mPairs and mIndices. This class works by preallocating
      memory, otherwise you would need to realloc memory every time one or more pairs are added.
      Alternatively, you could implement this class as a stack, that dynamically allocates memory.

      If you need to safe memory, allocate only as much as you need by giving the correct number
      of pairs to the constructor, or alternatively, create a big class and then create a copy 
      of the object with the copy constructor.

      When modifying the dots, the allocated size of the memory stays always the same.

 */

#define NODOT -1 

//------------------------------------< constructors and destructors >-----
ImplAlignmentMatrix::ImplAlignmentMatrix( long ndots ) : ImplAlignment(), 
mIndex(NULL),
mAllocatedIndexSize(0) 
{
	debug_func_cerr(5);

}

/*  It this important to add mPairs(NULL) to the list of initializer, because
	otherwise (at least gcc) adds an implicit mPairs(src.mPairs), which is 
	absolutely disastrous 
 */
ImplAlignmentMatrix::ImplAlignmentMatrix( const ImplAlignmentMatrix & src) : 
	ImplAlignment( src ), 
	mPairs(),			
	mIndex(NULL),
	mAllocatedIndexSize( src.mAllocatedIndexSize) 
	{
	debug_func_cerr(5);

	// create a deep copy of src.mPairs
	PairConstIterator it(src.mPairs.begin()), it_end(src.mPairs.end());
	for (; it != it_end; ++it) 
	{
		ResiduePAIR * p = *it;
		mPairs.push_back( new ResiduePAIR( *p ) );
	}

	// create a copy of the index (if existing)
	if (src.mIndex) 
	{
		mIndex = new Dot[mAllocatedIndexSize];
		memcpy( mIndex, src.mIndex, sizeof(Dot) * (mAllocatedIndexSize));
	}
	}

ImplAlignmentMatrix::~ImplAlignmentMatrix( ) 
{
	debug_func_cerr(5);

	clear();

	// mPairs and mIndex is deleted only now
	if (mIndex != NULL)
		delete [] mIndex;

	mIndex = NULL; 
}


//-----------------------------------------------------------------------------------------------------------   

AlignmentConstIterator ImplAlignmentMatrix::begin() const 
{ 
	if (mChangedLength) calculateLength();
	return AlignmentConstIterator( new ImplAlignmentMatrix_ConstIterator( mPairs, 0, mPairs.size() )); 
}

AlignmentConstIterator ImplAlignmentMatrix::end() const 
{ 
	if (mChangedLength) calculateLength();
	return AlignmentConstIterator( new ImplAlignmentMatrix_ConstIterator(mPairs, mPairs.size(), mPairs.size() )); 
}
//-----------------------------------------------------------------------------------------------------------   

AlignmentIterator ImplAlignmentMatrix::begin() 
{ 
	if (mChangedLength) calculateLength();
	return AlignmentIterator( new ImplAlignmentMatrix_Iterator( mPairs, 0, mPairs.size() )); 
}

AlignmentIterator ImplAlignmentMatrix::end() 
{ 
	if (mChangedLength) calculateLength();
	return AlignmentIterator( new ImplAlignmentMatrix_Iterator(mPairs, mPairs.size(), mPairs.size() )); 
}

//----------------> accessors <------------------------------------------------------------------------------

ResiduePAIR ImplAlignmentMatrix::front() const 
{ 
	if (mChangedLength) calculateLength(); 
	return (mPairs.size() > 0) ? (*mPairs.front()) : ResiduePAIR(NO_POS,NO_POS,0); 
}

ResiduePAIR ImplAlignmentMatrix::back()  const 
{ 
	if (mChangedLength) calculateLength(); 
	return (mPairs.size() > 0) ? (*mPairs.back()) : ResiduePAIR(NO_POS,NO_POS,0); 
}

//-------------------------------------------------------------------------------------------------------------
void ImplAlignmentMatrix::addPair( ResiduePAIR * pair ) 
{ 
	ImplAlignment::addPair( pair );
	mPairs.push_back( pair );
	setChangedLength();
} 

//-------------------------------------------------------------------------------------------------------------
/** retrieves a pair of residues from the alignment */
ResiduePAIR ImplAlignmentMatrix::getPair( const ResiduePAIR & p) const 
{
  /** generic implementation - returns any pair of row */
  PAIRVECTOR::iterator it(mPairs.begin()), it_end(mPairs.end());
  for (;it != it_end; ++it)
    if ((*it)->mRow == p.mRow)
      return **it;
  
  return ResiduePAIR();
} 

//-----------------------------------------------------------------------------------------------------------   
void ImplAlignmentMatrix::updateBoundaries() const
{

	Position max_size = mPairs.size();

	mRowFrom = mRowTo = mColFrom = mColTo = NO_POS;
	
	// ignore empty alignments
	if (max_size == 0)
	  return;
		
	mRowFrom = std::numeric_limits<Position>::max();
	mColFrom = std::numeric_limits<Position>::max();
	mRowTo = std::numeric_limits<Position>::min();
	mColTo = std::numeric_limits<Position>::min();
	
	PAIRVECTOR::const_iterator it(mPairs.begin()), it_end(mPairs.end());
	for (; it != it_end; ++it )
	{
    	const Position row = (*it)->mRow;
    	const Position col = (*it)->mCol;
    	
		// get maximum boundaries
    	if (row < mRowFrom) mRowFrom = row;
    	if (col < mColFrom) mColFrom = col;
    	if (row > mRowTo)   mRowTo = row;
    	if (col > mColTo)   mColTo = col;		
	}
	++mRowTo;
	++mColTo;
}

//----------------------------------------------------------------------------------------------------------------------
void ImplAlignmentMatrix::removePair( const ResiduePAIR & p ) 
{ 
	debug_func_cerr(5);

	PAIRVECTOR::iterator it(mPairs.begin());
	// re-check mPairs.end() at each loop iteration, do not use
	// test against other iterator, as this will not catch
	// the deletion of the last element.
	for ( ;it != mPairs.end(); )
	{
		if ( **it == p )
		{
			delete *it;
			it = mPairs.erase( it );
		}
		else
			++it;
	}

	updateBoundaries();
	setChangedLength();    
} 

//--------------------------------------------------------------------------------------------------------
void ImplAlignmentMatrix::clear() 
{
	debug_func_cerr(5);

	ImplAlignment::clear();

	mRowFrom = mRowTo = mColFrom = mColTo = NO_POS;
	if (mIndex != NULL)
		delete [] mIndex; 
	mIndex = NULL;

	PAIRVECTOR::iterator it(mPairs.begin()), it_end(mPairs.end());
	for (;it != it_end; ++it) delete *it;
	mPairs.clear();
}

//------------------------------------> sorting subroutines <-----------------------------------------------

//----------------------------------------------------------------------------------------------------------
void ImplAlignmentMatrix::sortDotsByDiagonal(Position from, Position to) const 
{
	Position lastsmall, c1, i;
	Diagonal m; /* value of median */
	ResiduePAIR * t;

	if (from < to) 
	{
		c1 = (from + to) / 2;

		t = mPairs[from]; mPairs[from] = mPairs[c1]; mPairs[c1] = t;

		m = calculateDiagonal( *mPairs[from] );					// choose median-value to compare to
		lastsmall = from;
		for (i = from + 1; i <= to; i++) 
		{
			if ( calculateDiagonal(*mPairs[i]) < m) 
			{						// swap lastsmall and i
				lastsmall++;
				t = mPairs[lastsmall]; mPairs[lastsmall] = mPairs[i]; mPairs[i] = t;
			}
		}

		t = mPairs[from]; mPairs[from] = mPairs[lastsmall]; mPairs[lastsmall] = t;

		sortDotsByDiagonal( from, lastsmall);
		sortDotsByDiagonal( lastsmall + 1, to);
	}
}

//----------------------------------------------------------------------------------------------------------------------
void ImplAlignmentMatrix::sortDotsByRow(Position from, Position to) const 
{
	Position lastsmall, c1, i;
	Position m; /* value of median */
	ResiduePAIR * t;

	if (from < to) 
	{
		c1 = (from + to) / 2;

		t = mPairs[from]; 
		mPairs[from] = mPairs[c1]; 
		mPairs[c1] = t;

		m = mPairs[from]->mRow;                                 // choose median-value to compare to
		lastsmall = from;
		for (i = from + 1; i <= to; i++) 
		{
			if ( mPairs[i]->mRow < m) 
			{                         // swap lastsmall and i
				lastsmall++;
				t = mPairs[lastsmall]; 
				mPairs[lastsmall] = mPairs[i]; 
				mPairs[i] = t;
			}
		}

		t = mPairs[from]; mPairs[from] = mPairs[lastsmall]; mPairs[lastsmall] = t;

		m = lastsmall;
		sortDotsByRow( from, m);
		sortDotsByRow( m+1, to);
	}
}

//----------------------------------------------------------------------------------------------------------------------
void ImplAlignmentMatrix::sortDotsByCol(Position from, Position to) const 
{
	Position lastsmall, c1, i;
	Position m; /* value of median */
	ResiduePAIR * t;

	if (from < to ) 
	{
		c1 = (from + to) / 2;

		t = mPairs[from]; 
		mPairs[from] = mPairs[c1]; 
		mPairs[c1] = t;

		m = mPairs[from]->mCol;                                 // choose median-value to compare to
		lastsmall = from;
		for (i = from + 1; i <= to; i++) 
		{
			if ( mPairs[i]->mCol < m) 
			{                         // swap lastsmall and i
				lastsmall++;
				t = mPairs[lastsmall]; 
				mPairs[lastsmall] = mPairs[i]; 
				mPairs[i] = t;
			}
		}

		t = mPairs[from]; mPairs[from] = mPairs[lastsmall]; mPairs[lastsmall] = t;

		m = lastsmall;
		sortDotsByCol( from, m);
		sortDotsByCol( m+1, to);
	}
}


//-------------------------------------------------------------------------------------------------------------------- 
void ImplAlignmentMatrix::eliminateDuplicates() const 
{

	mRowFrom = std::numeric_limits<Position>::max();
	mColFrom = std::numeric_limits<Position>::max();
	mRowTo = std::numeric_limits<Position>::min();
	mColTo = std::numeric_limits<Position>::min();

	Position last_row = NO_POS;
	Position last_col = NO_POS;

	std::vector<ResiduePAIR *> temp = mPairs;
	mPairs.clear();
	mPairs.reserve(temp.size());

	PairConstIterator it(temp.begin()), it_end(temp.end());

	for (; it != it_end; ++it) {

		Position row = (*it)->mRow;
		Position col = (*it)->mCol;

		/* delete pairs not needed any more */
		if (last_row == row && last_col == col) 
		{
			delete (*it);
			continue;
		}

		// get maximum boundaries
		if (row < mRowFrom) mRowFrom = row;
		if (col < mColFrom) mColFrom = col;
		if (row > mRowTo)   mRowTo = row;
		if (col > mColTo)   mColTo = col;

		mPairs.push_back(*it);
		last_row = row;
		last_col = col;
	}

	++mRowTo;
	++mColTo;
}


//--------------------------------------------------------------------------------------------------------------
void ImplAlignmentMatrix::allocateIndex( unsigned long size ) const 
{
	if (mIndex != NULL) 
		delete [] mIndex;

	mAllocatedIndexSize = size;
	mIndex = new Dot[size];
}

//--------------------------------------------------------------------------------------------------------------
// instead of calculating the length, it resorts the dots :=)
void ImplAlignmentMatrix::calculateLength() const 
{
	debug_func_cerr(5);

	mChangedLength = false;
	mRowFrom = mRowTo = mColFrom = mColTo = NO_POS;
	if (mIndex != NULL)
		delete [] mIndex;
	mIndex = NULL;

	mAllocatedIndexSize = 0;
	setNumGaps(0);
	setLength(mPairs.size());

	if (mPairs.empty()) 
		return;

	// 1. sort Dots
	sortDots();

	// 2. eliminate duplicates. At the same time, this sets mRowFrom, mRowTo, etc.
	eliminateDuplicates();

	// 3. build index for quick access (row, col, diagonal, etc.)
	buildIndex();

}  

} // namespace alignlib
