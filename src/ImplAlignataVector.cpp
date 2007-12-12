/*
  alignlib - a library for aligning protein sequences

  $Id: ImplAlignataVector.cpp,v 1.4 2004/06/02 12:11:37 aheger Exp $

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
#include "ImplAlignataVector.h"
#include "AlignataIterator.h"

#ifdef WITH_DMALLOC
#include <dmalloc.h>
#endif

using namespace std;

namespace alignlib {

/** by how much the vector grows */
#define GROWTH_FACTOR 2

//------------------------------factory functions -----------------------------
Alignata * makeAlignataVector() {
	return new ImplAlignataVector();
}

//------------------------------------< constructors and destructors >-----
ImplAlignataVector::ImplAlignataVector() : 
	ImplAlignata(), 
	mRowFrom(NO_POS), 
	mRowTo(NO_POS) {
}

ImplAlignataVector::ImplAlignataVector( const ImplAlignataVector& src) : 
	ImplAlignata( src ), 
	mRowFrom( src.mRowFrom), 
	mRowTo( src.mRowTo) {
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

ImplAlignataVector::~ImplAlignataVector( ) 
{
	debug_func_cerr(5);

	clear();
}

//------------------------------------------------------------------------------------------------------------
ImplAlignataVector * ImplAlignataVector::getNew() const {
	return new ImplAlignataVector();
}

ImplAlignataVector * ImplAlignataVector::getClone() const {
	return new ImplAlignataVector( *this );
}

//-----------------------------------------------------------------------------------------------------------   

AlignataConstIterator ImplAlignataVector::begin() const { 
	// skip gaps in the beginning of the alignment
	return AlignataConstIterator( new ImplAlignataVector_ConstIterator( mPairs, mRowFrom, mRowFrom, mRowTo )); 
}

AlignataConstIterator ImplAlignataVector::end() const { 
	return AlignataConstIterator( new ImplAlignataVector_ConstIterator(mPairs, NO_POS, mRowFrom, mRowTo)); 
}

//-----------------------------------------------------------------------------------------------------------   

AlignataIterator ImplAlignataVector::begin() { 
	// skip gaps in the beginning of the alignment
	return AlignataIterator( new ImplAlignataVector_Iterator( mPairs, mRowFrom, mRowFrom, mRowTo )); 
}

AlignataIterator ImplAlignataVector::end() { 
	return AlignataIterator( new ImplAlignataVector_Iterator(mPairs, NO_POS, mRowFrom, mRowTo)); 
}

//----------------> accessors <------------------------------------------------------------------------------

Position ImplAlignataVector::getRowFrom() const { return mRowFrom; }
Position ImplAlignataVector::getColFrom() const { 
	if (mRowFrom != NO_POS) 
		return mPairs[mRowFrom]->mCol; 
	else 
		return NO_POS;
}

Position ImplAlignataVector::getRowTo()   const { return mRowTo; }
Position ImplAlignataVector::getColTo()   const { 
	if (mRowTo != NO_POS) 
		return mPairs[mRowTo-1]->mCol+1; 
	else
		return NO_POS;
}


ResiduePAIR ImplAlignataVector::front() const { return *mPairs[mRowFrom]; }
ResiduePAIR ImplAlignataVector::back()  const { return *mPairs[mRowTo-1]; }

//----------------------------------------------------------------------------------------------------------------------
void ImplAlignataVector::addPair( ResiduePAIR * new_pair ) { 

	Position new_row = new_pair->mRow;

	if (mRowFrom > new_row || mRowFrom == NO_POS) 
		mRowFrom = new_row;

	if (mRowTo <= new_row) 
		mRowTo = new_row + 1;

	unsigned int needed_size = mRowTo;
	if (mPairs.size() < needed_size) 
		mPairs.resize( needed_size * GROWTH_FACTOR, NULL);

	if (mPairs[new_row] != NULL) 
		delete mPairs[new_row];

	mPairs[new_row] = new_pair; 
	setChangedLength();
} 

//--------------------------------------------------------------------------------------------------------------
void ImplAlignataVector::moveAlignment( Position row_offset, Position col_offset) 
{
	debug_func_cerr(5);


	// create copy of mPairs (= copy of pointers)
	PAIRVECTOR copy(mPairs);

	PairIterator it(copy.begin()), it_end(copy.end());

	// set new alignment coordinates
	mRowFrom += row_offset;  
	mRowTo   += row_offset;
	unsigned int needed_size = mRowTo;

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

}


//----------------------------------------------------------------------------------------------------------
/** retrieves a pair of residues from the alignment */
ResiduePAIR ImplAlignataVector::getPair( const ResiduePAIR & p) const 
{
	if (p.mRow != NO_POS)
		return *mPairs[p.mRow];
	else
		return ResiduePAIR();
} 

//----------------------------------------------------------------------------------------------------------
void ImplAlignataVector::removePair( const ResiduePAIR & old_pair ) { 

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
void ImplAlignataVector::clearContainer() { 

	PAIRVECTOR::iterator it(mPairs.begin()), it_end(mPairs.end());
	for (;it != it_end; ++it) delete *it;
	mPairs.clear();

}

//----------------------------------------------------------------------------------------------------------------
void ImplAlignataVector::clear() { 
	ImplAlignata::clear();

	clearContainer();

	mRowFrom = NO_POS;
	mRowTo = NO_POS; 
}

//--------------> mapping functions <----------------------------------------------------------------------------
Position ImplAlignataVector::mapRowToCol( Position pos, SearchType search ) const 
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
void ImplAlignataVector::resetBoundaries() {

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
void ImplAlignataVector::removeRowRegion( Position from, Position to) {

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
void ImplAlignataVector::removeColRegion( Position from, Position to) {

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

