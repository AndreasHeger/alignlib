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
ImplAlignataVector::ImplAlignataVector() : ImplAlignata(), mRowFrom(0), mRowTo(0) {
}

ImplAlignataVector::ImplAlignataVector( const ImplAlignataVector& src) : 
  ImplAlignata( src ), 
  mRowFrom( src.mRowFrom), 
  mRowTo( src.mRowTo) 
{
  debug_func_cerr(5);


    // do not call clear() function, as this will delete attributes of parent
    // class like mScore.

    clearContainer();
    mPairs.resize( mRowTo + 1, NULL);

    PAIRVECTOR::const_iterator it(src.mPairs.begin()), it_end(src.mPairs.end());
    for (;it != it_end; ++it) 
      if (*it != NULL) 
	mPairs[(*it)->mRow] = new ResiduePAIR( (**it) );
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
  return AlignataConstIterator( new ImplAlignataVector_ConstIterator(mPairs, 0, mRowFrom, mRowTo)); 
}

//-----------------------------------------------------------------------------------------------------------   

AlignataIterator ImplAlignataVector::begin() { 
  // skip gaps in the beginning of the alignment
  return AlignataIterator( new ImplAlignataVector_Iterator( mPairs, mRowFrom, mRowFrom, mRowTo )); 
}

AlignataIterator ImplAlignataVector::end() { 
  return AlignataIterator( new ImplAlignataVector_Iterator(mPairs, 0, mRowFrom, mRowTo)); 
}

//----------------> accessors <------------------------------------------------------------------------------

Position ImplAlignataVector::getRowFrom() const { return mRowFrom; }
Position ImplAlignataVector::getColFrom() const { 
  if (mRowFrom) 
    return mPairs[mRowFrom]->mCol; 
  else 
    return 0;
}

Position ImplAlignataVector::getRowTo()   const { return mRowTo; }
Position ImplAlignataVector::getColTo()   const { 
  if (mRowTo != 0) 
    return mPairs[mRowTo]->mCol; 
  else
    return 0;
}

  
ResiduePAIR ImplAlignataVector::getFirstPair() const { return *mPairs[mRowFrom]; }
ResiduePAIR ImplAlignataVector::getLastPair()  const { return *mPairs[mRowTo]; }

//----------------------------------------------------------------------------------------------------------------------
void ImplAlignataVector::addPair( ResiduePAIR * new_pair ) { 
  
  Position new_row = new_pair->mRow;
  

  if (mRowFrom > new_row || mRowFrom == 0) 
    mRowFrom = new_row;

  if (mRowTo < new_row) 
    mRowTo = new_row;

  unsigned int needed_size = mRowTo + 1;
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
  unsigned int needed_size = mRowTo + 1;

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
ResiduePAIR ImplAlignataVector::getPair( Position row ) const {
  return *mPairs[row];
} 

//----------------------------------------------------------------------------------------------------------
void ImplAlignataVector::removePair( const ResiduePAIR & old_pair ) { 

  // no resizing is done;
  if (old_pair.mRow >= mRowFrom && old_pair.mRow <= mRowTo)
    mPairs[old_pair.mRow] = NULL;

  // if first pair has been deleted, update mRowFrom
  while (mRowFrom <= mRowTo && mPairs[mRowFrom] == NULL) 
    mRowFrom++;
  
  // if last pairs has been deleted, update mRowTo
  while (mRowTo >= mRowFrom && mPairs[mRowTo] == NULL)
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
  
  mRowFrom = 0;
  mRowTo = 0; 
}

//--------------> mapping functions <----------------------------------------------------------------------------
Position ImplAlignataVector::mapRowToCol( Position pos ) const {
  if (pos >= mRowFrom && pos <= mRowTo && mPairs[pos] != NULL)
    return mPairs[pos]->mCol;
  else
    return 0;
}


//-----------------------------------------------------------------------------------------------------------   
void ImplAlignataVector::resetBoundaries() {

  Position max_size = mPairs.size();

  // find new mRowFrom and mRowTo, if they have changed
  while (mRowFrom < max_size && mPairs[mRowFrom] == NULL) mRowFrom++;
  while (mRowTo > 0          && mPairs[mRowTo] == NULL) mRowTo--;

  // test >= because max_size might be zero, if all has been 
  // deleted.
  if (mRowFrom >= max_size) 
    mRowFrom = 0;
}
//-----------------------------------------------------------------------------------------------------------   
void ImplAlignataVector::removeRowRegion( Position from, Position to) {

  Position pos;

  if (from < mRowFrom)
    from = mRowFrom;
  
  if (to > mRowTo)
    to = mRowTo;
  
  // delete aligned positions
  for ( pos = from; pos <= to; pos++) {
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

  for (pos = mRowFrom; pos <= mRowTo; pos++)
    if (mPairs[pos] != NULL && mPairs[pos]->mCol >= from && mPairs[pos]->mCol <= to) {
      delete mPairs[pos];
      mPairs[pos] = NULL;
    }

  resetBoundaries();

  setChangedLength();

  return;
}




} // namespace alignlib

