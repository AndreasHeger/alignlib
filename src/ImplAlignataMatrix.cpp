/*
  alignlib - a library for aligning protein sequences

  $Id: ImplAlignataMatrix.cpp,v 1.6 2004/09/24 19:03:27 aheger Exp $

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
#include "ImplAlignataMatrix.h"
#include "AlignataIterator.h"
#include "AlignException.h"

#ifdef WITH_DMALLOC
#include <dmalloc.h>
#endif

using namespace std;

namespace alignlib {

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
  ImplAlignataMatrix::ImplAlignataMatrix( long ndots ) : ImplAlignata(), 
  mIndex(NULL),
  mRowFrom(NO_POS), mRowTo(NO_POS), mColFrom(NO_POS), mColTo(NO_POS), 
  mAllocatedIndexSize(0) 
  {
    debug_func_cerr(5);

  }

  /*  It this important to add mPairs(NULL) to the list of initializer, because
	otherwise (at least gcc) adds an implicit mPairs(src.mPairs), which is 
	absolutely disastrous 
   */
  ImplAlignataMatrix::ImplAlignataMatrix( const ImplAlignataMatrix & src) : 
    ImplAlignata( src ), 
    mPairs(),			
    mIndex(NULL),
    mRowFrom( src.mRowFrom), 
    mRowTo( src.mRowTo),
    mColFrom( src.mColFrom), 
    mColTo( src.mColTo),
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

  ImplAlignataMatrix::~ImplAlignataMatrix( ) 
    {
      debug_func_cerr(5);

      clear();

      // mPairs and mIndex is deleted only now
      if (mIndex != NULL)
        delete [] mIndex;

      mIndex = NULL; 
    }


  //-----------------------------------------------------------------------------------------------------------   

  AlignataConstIterator ImplAlignataMatrix::begin() const { 
    if (mChangedLength) calculateLength();
    return AlignataConstIterator( new ImplAlignataMatrix_ConstIterator( mPairs, 0, mPairs.size() )); 
  }

  AlignataConstIterator ImplAlignataMatrix::end() const { 
    if (mChangedLength) calculateLength();
    return AlignataConstIterator( new ImplAlignataMatrix_ConstIterator(mPairs, mPairs.size(), mPairs.size() )); 
  }
  //-----------------------------------------------------------------------------------------------------------   

  AlignataIterator ImplAlignataMatrix::begin() { 
    if (mChangedLength) calculateLength();
    return AlignataIterator( new ImplAlignataMatrix_Iterator( mPairs, 0, mPairs.size() )); 
  }

  AlignataIterator ImplAlignataMatrix::end() { 
    if (mChangedLength) calculateLength();
    return AlignataIterator( new ImplAlignataMatrix_Iterator(mPairs, mPairs.size(), mPairs.size() )); 
  }

  //----------------> accessors <------------------------------------------------------------------------------

  Position ImplAlignataMatrix::getRowFrom() const { if (mChangedLength) calculateLength(); return mRowFrom; }
  Position ImplAlignataMatrix::getColFrom() const { if (mChangedLength) calculateLength(); return mColFrom; }
  Position ImplAlignataMatrix::getRowTo()   const { if (mChangedLength) calculateLength(); return mRowTo; }
  Position ImplAlignataMatrix::getColTo()   const { if (mChangedLength) calculateLength(); return mColTo; }

  ResiduePAIR ImplAlignataMatrix::front() const 
  { 
    if (mChangedLength) calculateLength(); 
    return (mPairs.size() > 0) ? (*mPairs.front()) : ResiduePAIR(NO_POS,NO_POS,0); 
  }

  ResiduePAIR ImplAlignataMatrix::back()  const 
  { 
    if (mChangedLength) calculateLength(); 
    return (mPairs.size() > 0) ? (*mPairs.back()) : ResiduePAIR(NO_POS,NO_POS,0); 
  }

  //-------------------------------------------------------------------------------------------------------------
  void ImplAlignataMatrix::addPair( ResiduePAIR * new_pair ) 
  { 

    mPairs.push_back( new_pair );
    setChangedLength();		// has to be resorted
  } 
  //-------------------------------------------------------------------------------------------------------------
  void ImplAlignataMatrix::addPair( Position row, Position col, Score score ) 
  { 

    mPairs.push_back( new ResiduePAIR(row,col,score) );
    setChangedLength();		// has to be resorted
  } 

  //-------------------------------------------------------------------------------------------------------------
  /** retrieves a pair of residues from the alignment */
  ResiduePAIR ImplAlignataMatrix::getPair( Position row ) const 
  {
    // TODO!! to be implemented
	assert( false );
    return *mPairs[row];
  } 

  //----------------------------------------------------------------------------------------------------------------------
  void ImplAlignataMatrix::removePair( const ResiduePAIR & old_pair ) 
  { 

    setChangedLength();
    assert( false );
    // TODO !! to be implemented.  move and split highly inefficient

  } 

  //--------------------------------------------------------------------------------------------------------
  void ImplAlignataMatrix::clear() 
    {
      debug_func_cerr(5);

      ImplAlignata::clear();

      mRowFrom = NO_POS;
      mRowTo = NO_POS;
      mColFrom = NO_POS;
      mRowTo = NO_POS;
      if (mIndex != NULL)
        delete [] mIndex; 
      mIndex = NULL;

      PAIRVECTOR::iterator it(mPairs.begin()), it_end(mPairs.end());
      for (;it != it_end; ++it) delete *it;
      mPairs.clear();
    }

  //--------------------------------------------------------------------------------------------------------------
  void ImplAlignataMatrix::moveAlignment( Position row_offset, Position col_offset) 
    {
      debug_func_cerr(5);

      ImplAlignata::moveAlignment( row_offset, col_offset);

      // reset coordinates of alignment
      mRowFrom += row_offset;
      mRowTo   += row_offset;
      mColFrom += col_offset;
      mColTo   += col_offset;

    }
  //------------------------------------> sorting subroutines <-----------------------------------------------

  //----------------------------------------------------------------------------------------------------------
  void ImplAlignataMatrix::sortDotsByDiagonal(Position from, Position to) const 
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
  void ImplAlignataMatrix::sortDotsByRow(Position from, Position to) const 
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
  void ImplAlignataMatrix::sortDotsByCol(Position from, Position to) const 
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
      for (i = from + 1; i <= to; i++) {
        if ( mPairs[i]->mCol < m) {                         // swap lastsmall and i
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
  void ImplAlignataMatrix::eliminateDuplicates() const {

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
      if (last_row == row && last_col == col) {
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
  void ImplAlignataMatrix::allocateIndex( unsigned long size ) const 
  {
    if (mIndex != NULL) 
      delete [] mIndex;

    mAllocatedIndexSize = size;
    mIndex = new Dot[size];
  }

  //--------------------------------------------------------------------------------------------------------------
  // instead of calculating the length, it resorts the dots :=)
  void ImplAlignataMatrix::calculateLength() const 
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
