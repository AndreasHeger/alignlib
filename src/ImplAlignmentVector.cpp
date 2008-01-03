/*
  alignlib - a library for aligning protein sequences

  $Id: ImplAlignmentVector.cpp,v 1.4 2004/06/02 12:11:37 aheger Exp $

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
#include "alignlib.h"
#include "AlignlibDebug.h"
#include "ImplAlignmentVector.h"
#include "AlignmentIterator.h"

using namespace std;

namespace alignlib 
{

/** by how much the vector grows */
#define GROWTH_FACTOR 2

//------------------------------factory functions -----------------------------
Alignment * makeAlignmentVector() 
{
	return new ImplAlignmentVector();
}

//------------------------------------< constructors and destructors >-----
ImplAlignmentVector::ImplAlignmentVector() : 
	ImplAlignment()
{
}

ImplAlignmentVector::ImplAlignmentVector( const ImplAlignmentVector& src) : 
	ImplAlignment( src )
{
	debug_func_cerr(5);

	// do not call clear() function, as this will delete attributes of parent
	// class like mScore.

	clearContainer();
	if (mRowTo == NO_POS)
		mPairs.clear();
	else
	{
		mPairs.resize( mRowTo, NULL);

		PAIRVECTOR::const_iterator it(src.mPairs.begin()), it_end(src.mPairs.end());
		for (;it != it_end; ++it) 
			if (*it != NULL) 
				mPairs[(*it)->mRow] = new ResiduePAIR( (**it) );
	}
}

ImplAlignmentVector::~ImplAlignmentVector( ) 
{
	debug_func_cerr(5);

	clear();
}

//------------------------------------------------------------------------------------------------------------
ImplAlignmentVector * ImplAlignmentVector::getNew() const
{
	return new ImplAlignmentVector();
}

ImplAlignmentVector * ImplAlignmentVector::getClone() const 
{
	return new ImplAlignmentVector( *this );
}

//-----------------------------------------------------------------------------------------------------------   

AlignmentConstIterator ImplAlignmentVector::begin() const 
{ 	
	return AlignmentConstIterator( new ImplAlignmentVector_ConstIterator( mPairs, mRowFrom, mRowFrom, mRowTo )); 
}

AlignmentConstIterator ImplAlignmentVector::end() const 
{ 
	return AlignmentConstIterator( new ImplAlignmentVector_ConstIterator(mPairs, NO_POS, mRowFrom, mRowTo)); 
}

//-----------------------------------------------------------------------------------------------------------   

AlignmentIterator ImplAlignmentVector::begin() 
{	
	return AlignmentIterator( new ImplAlignmentVector_Iterator( mPairs, mRowFrom, mRowFrom, mRowTo )); 
}

AlignmentIterator ImplAlignmentVector::end() 
{ 
	return AlignmentIterator( new ImplAlignmentVector_Iterator(mPairs, NO_POS, mRowFrom, mRowTo)); 
}

//----------------> accessors <------------------------------------------------------------------------------

ResiduePAIR ImplAlignmentVector::front() const { return *mPairs[mRowFrom]; }
ResiduePAIR ImplAlignmentVector::back()  const { return *mPairs[mRowTo-1]; }

//----------------------------------------------------------------------------------------------------------------------
void ImplAlignmentVector::addPair( ResiduePAIR * new_pair ) 
{ 

	ImplAlignment::addPair( new_pair );
	
	Position new_row = new_pair->mRow;

	size_t needed_size = std::max( mPairs.size(), (size_t)new_row + 1);
	
	if (mPairs.size() < needed_size) 
		mPairs.resize( needed_size * GROWTH_FACTOR, NULL);

	if (mPairs[new_row] != NULL) 
		delete mPairs[new_row];

	mPairs[new_row] = new_pair; 
	setChangedLength();
} 

//--------------------------------------------------------------------------------------------------------------
void ImplAlignmentVector::moveAlignment( Position row_offset, Position col_offset) 
{
	debug_func_cerr(5);

	// create copy of mPairs (= copy of pointers)
	PAIRVECTOR copy(mPairs);

	PairIterator it(copy.begin()), it_end(copy.end());

	size_t needed_size = std::max( mPairs.size(), (size_t)mRowTo + row_offset );

	// delete old alignment and allocate needed size
	mPairs.clear();
	mPairs.resize( needed_size, NULL);

	// copy pointers from copy into mPairs
	for (; it != it_end; ++it) {
		ResiduePAIR * p = *it;
		if (p != NULL) {
			p->mRow += row_offset;
			p->mCol += col_offset;
			mPairs[p->mRow] = p;
		}
	}
	// set new alignment coordinates
	mRowFrom += row_offset;  
	mRowTo   += row_offset;
	mColFrom += col_offset;
	mColTo += col_offset;

}


//----------------------------------------------------------------------------------------------------------
/** retrieves a pair of residues from the alignment */
ResiduePAIR ImplAlignmentVector::getPair( const ResiduePAIR & p) const 
{
	if (p.mRow != NO_POS)
		return *mPairs[p.mRow];
	else
		return ResiduePAIR();
} 

//----------------------------------------------------------------------------------------------------------
void ImplAlignmentVector::removePair( const ResiduePAIR & old_pair ) 
{ 

	// no resizing is done;
	if (old_pair.mRow >= mRowFrom && old_pair.mRow < mRowTo)
		mPairs[old_pair.mRow] = NULL;

	// if first pair has been deleted, update mRowFrom
	while (mRowFrom < mRowTo && mPairs[mRowFrom] == NULL) 
		mRowFrom++;

	// if last pairs has been deleted, update mRowTo
	while (mRowTo >= mRowFrom && mPairs[mRowTo-1] == NULL)
		mRowTo--;

	if (mRowFrom > mRowTo)
		clear();

	setChangedLength();
} 



//----------------------------------------------------------------------------------------------------------------
void ImplAlignmentVector::clearContainer() 
{ 

	PAIRVECTOR::iterator it(mPairs.begin()), it_end(mPairs.end());
	for (;it != it_end; ++it) delete *it;
	mPairs.clear();

}

//----------------------------------------------------------------------------------------------------------------
void ImplAlignmentVector::clear() 
{ 
	ImplAlignment::clear();

	clearContainer();
}

//--------------> mapping functions <----------------------------------------------------------------------------
Position ImplAlignmentVector::mapRowToCol( Position pos, SearchType search ) const 
{
  if (mRowFrom == NO_POS) return NO_POS;
  
  if ( search == LEFT && pos >= mRowTo)
    return mPairs[mRowTo - 1]->mCol;

  if ( search == RIGHT && pos < mRowFrom )
    return mPairs[mRowFrom]->mCol;

  if (pos < mRowFrom || pos >= mRowTo) return NO_POS;

	if (mPairs[pos] != NULL) 
		return mPairs[pos]->mCol;

	if (search == NO_SEARCH)
	{
		return NO_POS;
	}
	else if (search == LEFT)
	{
		while (pos >= mRowFrom && mPairs[pos] == NULL )
			--pos;
		if (pos < mRowFrom)
			return NO_POS;
	} 
	else if (search == RIGHT)
	{
		while (pos < mRowTo && mPairs[pos] == NULL )
			++pos;
		if (pos >= mRowTo)
			return NO_POS;
	}

	return mPairs[pos]->mCol;
}


//-----------------------------------------------------------------------------------------------------------   
void ImplAlignmentVector::resetBoundaries() {

	Position max_size = mPairs.size();

	// ignore empty alignments
	if (mRowFrom == NO_POS || max_size == 0)
	  return;
	
	// find new mRowFrom and mRowTo, if they have changed	
	while (mRowFrom < max_size && mPairs[mRowFrom] == NULL) ++mRowFrom;
	
        // test >= because max_size might be zero, if all has been 
        // deleted.
        if (mRowFrom >= max_size)
          {
            mRowFrom = NO_POS;
            mRowTo == NO_POS;
            return;
          }
	
	Position x = mRowTo - 1;
	while (x >= 0 && mPairs[x] == NULL) --x;
	mRowTo = x + 1;
	
}
//-----------------------------------------------------------------------------------------------------------   
void ImplAlignmentVector::removeRowRegion( Position from, Position to) {

	Position pos;

	if (from < mRowFrom)
		from = mRowFrom;

	if (to > mRowTo)
		to = mRowTo;

	// delete aligned positions
	for ( pos = from; pos < to; pos++) {
		if (mPairs[pos] != NULL) {
			delete mPairs[pos];
			mPairs[pos] = NULL;
		}
	}

	resetBoundaries();

	setChangedLength();

	return;
}

//-----------------------------------------------------------------------------------------------------------   
/* It is necessary to iterate from mFrowFrom to mRowTo, since the alignment need not be linear
 */
void ImplAlignmentVector::removeColRegion( Position from, Position to) {

	Position pos; 

	if (mRowFrom == NO_POS) return;
	
	for (pos = mRowFrom; pos < mRowTo; pos++)
		if (mPairs[pos] != NULL && mPairs[pos]->mCol >= from && mPairs[pos]->mCol < to) {
			delete mPairs[pos];
			mPairs[pos] = NULL;
		}

	resetBoundaries();

	setChangedLength();

	return;
}




} // namespace alignlib

