/*
  alignlib - a library for aligning protein sequences

  $Id: ImplAlignataHashDiagonal.cpp,v 1.4 2004/06/02 12:11:37 aheger Exp $

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
#include "ImplAlignataHashDiagonal.h"
#include "AlignataIterator.h"

#ifdef WITH_DMALLOC
#include <dmalloc.h>
#endif

using namespace std;

namespace alignlib {

  //------------------------------factory functions -----------------------------
  Alignata * makeAlignataHashDiagonal() 
  {
    return new ImplAlignataHashDiagonal();
  }

  //------------------------------------< constructors and destructors >-----
  ImplAlignataHashDiagonal::ImplAlignataHashDiagonal() : ImplAlignata() 
  {
    debug_func_cerr(5);

  }

  ImplAlignataHashDiagonal::ImplAlignataHashDiagonal( const ImplAlignataHashDiagonal& src) : ImplAlignata( src ) 
  {
    debug_func_cerr(5);


    clearContainer();

    PairIterator it(src.mPairs.begin()), it_end(src.mPairs.end());

    for (;it != it_end; ++it) 
      mPairs.insert( new ResiduePAIR(**it) );
  }

  ImplAlignataHashDiagonal::~ImplAlignataHashDiagonal( ) 
    {
      debug_func_cerr(5);

      clear();
    }

  //------------------------------------------------------------------------------------------------------------
  ImplAlignataHashDiagonal * ImplAlignataHashDiagonal::getNew() const 
  {
    return new ImplAlignataHashDiagonal();
  }

  ImplAlignataHashDiagonal * ImplAlignataHashDiagonal::getClone() const 
  {
    return new ImplAlignataHashDiagonal( *this );
  }

  //-------------------------------------------------------------------
  void ImplAlignataHashDiagonal::clearContainer() 
    { 
      PAIRSET::iterator it(mPairs.begin()), it_end(mPairs.end());
      for (;it != it_end; ++it) delete *it;
      mPairs.clear(); 
    }

  //-------------------------------------------------------------------
  void ImplAlignataHashDiagonal::clear() { 
    ImplAlignata::clear();

    clearContainer();
  }

  //-----------------------------------------------------------------------------------------------------------   

  AlignataConstIterator ImplAlignataHashDiagonal::begin() const 
  { 
    return AlignataConstIterator( new ImplAlignataHashDiagonal_ConstIterator(mPairs.begin()));
  }

  AlignataConstIterator ImplAlignataHashDiagonal::end()   const 
  { 
    return AlignataConstIterator( new ImplAlignataHashDiagonal_ConstIterator(mPairs.end())); 
  }

  AlignataIterator ImplAlignataHashDiagonal::begin() 
    { 
      return AlignataIterator( new ImplAlignataHashDiagonal_Iterator(mPairs.begin()));
    }

  AlignataIterator ImplAlignataHashDiagonal::end() 
    { 
      return AlignataIterator( new ImplAlignataHashDiagonal_Iterator(mPairs.end())); 
    }


  //----------------> accessors <------------------------------------------------------------------------------

  Position ImplAlignataHashDiagonal::getRowFrom() const { if (mChangedLength) calculateLength(); return mRowFrom; }
  Position ImplAlignataHashDiagonal::getColFrom() const { if (mChangedLength) calculateLength(); return mColFrom; }
  Position ImplAlignataHashDiagonal::getRowTo()   const { if (mChangedLength) calculateLength(); return mRowTo; }
  Position ImplAlignataHashDiagonal::getColTo()   const { if (mChangedLength) calculateLength(); return mColTo; }

  ResiduePAIR ImplAlignataHashDiagonal::front() const { return **(mPairs.begin()); }
  ResiduePAIR ImplAlignataHashDiagonal::back()  const { return **(mPairs.rbegin()); }

  void ImplAlignataHashDiagonal::addPair( ResiduePAIR * new_pair ) 
    {
      debug_func_cerr(5);

      if (mPairs.find( new_pair) != mPairs.end()) 
        {
          delete new_pair;
        } else 
          {
            setChangedLength(); 
            mPairs.insert( new_pair ); 
          } 
    } 

  //----------------------------------------------------------------------------------------------------------
  /** retrieves a pair of residues from the alignment */
  ResiduePAIR ImplAlignataHashDiagonal::getPair( const ResiduePAIR & p) const 
  {
	PairIterator it(mPairs.find( &p ));
	if (p != mPairs.end())
		return **it;
	else
		return ResiduePAIR();
  } 

  void ImplAlignataHashDiagonal::removePair( const ResiduePAIR & old_pair ) 
    {
      debug_func_cerr(5);

      setChangedLength(); 
      ResiduePAIR p(old_pair);
      mPairs.erase(&p); 
    } 

  void ImplAlignataHashDiagonal::calculateLength() const {

    mChangedLength = false;

    mRowFrom = std::numeric_limits<Position>::max();
    mRowTo = std::numeric_limits<Position>::min();
    mColFrom = std::numeric_limits<Position>::max();
    mColTo = std::numeric_limits<Position>::min();

    PairConstIterator it(mPairs.begin()), it_end(mPairs.end());

    Position length = 0;

    for (; it != it_end; ++it) 
      {
        Position col = (*it)->mCol;
        Position row = (*it)->mRow;

        // get maximum boundaries
        if (row < mRowFrom) mRowFrom = row;
        if (col < mColFrom) mColFrom = col;
        if (row > mRowTo)   mRowTo = row;
        if (col > mColTo)   mColTo = col;
        length++;
      }

    setNumGaps(0);
    setLength( length );
    ++mRowTo;
    ++mColTo;
  }

} // namespace alignlib
