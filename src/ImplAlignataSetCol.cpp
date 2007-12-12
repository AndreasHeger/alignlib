/*
  alignlib - a library for aligning protein sequences

  $Id: ImplAlignataSetCol.cpp,v 1.4 2004/06/02 12:11:37 aheger Exp $

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
#include "ImplAlignataSetCol.h"
#include "AlignataIterator.h"

#ifdef WITH_DMALLOC
#include <dmalloc.h>
#endif

using namespace std;

namespace alignlib {

  //------------------------------factory functions -----------------------------
  Alignata * makeAlignataSetCol() {
    return new ImplAlignataSetCol();
  }

  //------------------------------------< constructors and destructors >-----
  ImplAlignataSetCol::ImplAlignataSetCol() : ImplAlignata() 
  {
    debug_func_cerr(5);

  }

  ImplAlignataSetCol::ImplAlignataSetCol( const ImplAlignataSetCol& src) : ImplAlignata( src ) 
  {
    debug_func_cerr(5);


    clearContainer();

    PairConstIterator it(src.mPairs.begin()), it_end(src.mPairs.end());
    for (; it != it_end; ++it) 
      mPairs.insert( new ResiduePAIR(**it) );
  }

  ImplAlignataSetCol::~ImplAlignataSetCol( ) 
    {
      debug_func_cerr(5);

      clear();
    }

  //------------------------------------------------------------------------------------------------------------
  ImplAlignataSetCol * ImplAlignataSetCol::getNew() const {
    return new ImplAlignataSetCol();
  }

  ImplAlignataSetCol * ImplAlignataSetCol::getClone() const {
    return new ImplAlignataSetCol( *this );
  }

  //-------------------------------------------------------------------
  void ImplAlignataSetCol::clearContainer() { 
    PAIRSET::iterator it(mPairs.begin()), it_end(mPairs.end());
    for (;it != it_end; ++it) delete *it;
    mPairs.clear(); 
  }

  //------------------------------------------------------------------------------------------------------------
  void ImplAlignataSetCol::clear() { 
    ImplAlignata::clear();
    clearContainer();
  }

  //-----------------------------------------------------------------------------------------------------------   

  AlignataConstIterator ImplAlignataSetCol::begin() const { 
    return AlignataConstIterator( new ImplAlignataSetCol_ConstIterator(mPairs.begin()));
  }

  AlignataConstIterator ImplAlignataSetCol::end()   const { 
    return AlignataConstIterator( new ImplAlignataSetCol_ConstIterator(mPairs.end())); 
  }

  AlignataIterator ImplAlignataSetCol::begin() { 
    return AlignataIterator( new ImplAlignataSetCol_Iterator(mPairs.begin()));
  }

  AlignataIterator ImplAlignataSetCol::end() { 
    return AlignataIterator( new ImplAlignataSetCol_Iterator(mPairs.end())); 
  }


  //----------------> accessors <------------------------------------------------------------------------------

  Position ImplAlignataSetCol::getRowFrom() const { return (mPairs.size() > 0) ? (*mPairs.begin())->mRow : NO_POS; }
  Position ImplAlignataSetCol::getColFrom() const { return (mPairs.size() > 0) ? (*mPairs.begin())->mCol : NO_POS; }
  Position ImplAlignataSetCol::getRowTo()   const { return (mPairs.size() > 0) ? (*mPairs.rbegin())->mRow : NO_POS; }
  Position ImplAlignataSetCol::getColTo()   const { return (mPairs.size() > 0) ? (*mPairs.rbegin())->mCol : NO_POS; }

  ResiduePAIR ImplAlignataSetCol::front() const { return **(mPairs.begin()); }
  ResiduePAIR ImplAlignataSetCol::back()  const { return **(mPairs.rbegin()); }

  void ImplAlignataSetCol::addPair( ResiduePAIR * new_pair ) 
    {
      debug_func_cerr(5);

      if (mPairs.find( new_pair) != mPairs.end()) {
        delete new_pair;
      } else {
        setChangedLength(); 
        mPairs.insert( new_pair ); 
      } 
    } 

  //----------------------------------------------------------------------------------------------------------
  /** retrieves a pair of residues from the alignment */
  ResiduePAIR ImplAlignataSetCol::getPair( const ResiduePAIR & p) const {

    ResiduePAIR p(row, NO_POS, 0);
    PairIterator it(mPairs.find( &p ));
    return **it;
  } 


  void ImplAlignataSetCol::removePair( const ResiduePAIR & old_pair ) 
    {
      debug_func_cerr(5);


      setChangedLength(); 
      ResiduePAIR p(old_pair);
      mPairs.erase(&p); 
    } 

  //----------------------------------------------------------------------------------------
  /** This is non-generic routine. Since pairs are accessed by row, this is quite quick.
   */
  void ImplAlignataSetCol::removeColRegion( Position from, Position to) 
    {

      for (Position pos = from; pos < to; pos++) 
        {
          ResiduePAIR p(NO_POS, pos, 0);
          PairIterator it(mPairs.find( &p ));

          if (it != mPairs.end()) {
            setChangedLength(); 
            delete *it;
            mPairs.erase(it);
          }
        }

    }

  //----------------------------------------------------------------------------------------
  /** This is a generic routine. It creates a new alignment by making a copy of the old one */
  void ImplAlignataSetCol::removeRowRegion( Position from, Position to) 
    {

      PairIterator it(mPairs.begin()), it_end(mPairs.end());

      while (it != it_end)
        {
          if ( (*it)->mRow >= from && (*it)->mRow < to) 
            {
              PairIterator it2(it);
              ++it;
              delete *it2;
              mPairs.erase(it2);
            }
          else
            ++it;
        }

      setChangedLength();

    }


} // namespace alignlib

































