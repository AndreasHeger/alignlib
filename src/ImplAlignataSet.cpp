/*
  alignlib - a library for aligning protein sequences

  $Id: ImplAlignataSet.cpp,v 1.4 2004/06/02 12:11:37 aheger Exp $

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
#include "ImplAlignataSet.h"
#include "AlignataIterator.h"

#ifdef WITH_DMALLOC
#include <dmalloc.h>
#endif

using namespace std;

namespace alignlib {

  //------------------------------factory functions -----------------------------
  Alignata * makeAlignataSet() {
    return new ImplAlignataSet();
  }

  //------------------------------------< constructors and destructors >-----
  ImplAlignataSet::ImplAlignataSet() : ImplAlignata() 
  {
    debug_func_cerr(5);

  }

  ImplAlignataSet::ImplAlignataSet( const ImplAlignataSet& src) : ImplAlignata( src ) 
  {
    debug_func_cerr(5);


    clearContainer();

    PairIterator it(src.mPairs.begin()), it_end(src.mPairs.end());
    for (;it != it_end; ++it) 
      mPairs.insert( new ResiduePAIR(**it) );
  }

  ImplAlignataSet::~ImplAlignataSet( ) 
    {
      debug_func_cerr(5);

      clear();
    }

  //------------------------------------------------------------------------------------------------------------
  ImplAlignataSet * ImplAlignataSet::getNew() const {
    return new ImplAlignataSet();
  }

  ImplAlignataSet * ImplAlignataSet::getClone() const {
    return new ImplAlignataSet( *this );
  }

  //------------------------------------------------------------------------------------------------------------
  void ImplAlignataSet::clearContainer() { 
    PAIRSET::iterator it(mPairs.begin()), it_end(mPairs.end());
    for (;it != it_end; ++it) delete *it;
    mPairs.clear(); 
  }

  //------------------------------------------------------------------------------------------------------------
  void ImplAlignataSet::clear() { 
    ImplAlignata::clear();
    clearContainer();
  }

  //-----------------------------------------------------------------------------------------------------------   

  AlignataConstIterator ImplAlignataSet::begin() const { 
    return AlignataConstIterator( new ImplAlignataSet_ConstIterator(mPairs.begin()));
  }

  AlignataConstIterator ImplAlignataSet::end()   const { 
    return AlignataConstIterator( new ImplAlignataSet_ConstIterator(mPairs.end())); 
  }

  AlignataIterator ImplAlignataSet::begin() { 
    return AlignataIterator( new ImplAlignataSet_Iterator(mPairs.begin()));
  }

  AlignataIterator ImplAlignataSet::end() { 
    return AlignataIterator( new ImplAlignataSet_Iterator(mPairs.end())); 
  }


  //----------------> accessors <------------------------------------------------------------------------------

  Position ImplAlignataSet::getRowFrom() const { return (mPairs.size() > 0) ? (*mPairs.begin())->mRow : NO_POS; }
  Position ImplAlignataSet::getColFrom() const { return (mPairs.size() > 0) ? (*mPairs.begin())->mCol : NO_POS; }
  Position ImplAlignataSet::getRowTo()   const { return (mPairs.size() > 0) ? (*mPairs.rbegin())->mRow : NO_POS; }
  Position ImplAlignataSet::getColTo()   const { return (mPairs.size() > 0) ? (*mPairs.rbegin())->mCol : NO_POS; }

  ResiduePAIR ImplAlignataSet::front() const { return **(mPairs.begin()); }
  ResiduePAIR ImplAlignataSet::back()  const { return **(mPairs.rbegin()); }

  void ImplAlignataSet::addPair( ResiduePAIR * new_pair ) 
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
  ResiduePAIR ImplAlignataSet::getPair( Position row ) const {

    ResiduePAIR p(row, NO_POS, 0);
    PairIterator it(mPairs.find( &p ));
    return **it;
  } 


  void ImplAlignataSet::removePair( const ResiduePAIR & old_pair ) 
    {
      debug_func_cerr(5);


      ResiduePAIR p(old_pair);
      PairIterator it(mPairs.find( &p ));

      if (it != mPairs.end()) {
        setChangedLength(); 
        delete *it;
        mPairs.erase(it);
      }
    } 

  //----------------------------------------------------------------------------------------
  /** This is non-generic routine. Since pairs are accessed by row, this is quite quick.
   */
  void ImplAlignataSet::removeRowRegion( Position from, Position to) {

    for (Position pos = from; pos < to; pos++) {
      ResiduePAIR p(pos, NO_POS, 0);
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
  void ImplAlignataSet::removeColRegion( Position from, Position to) 
    {

      PairIterator it(mPairs.begin()), it_end(mPairs.end());

      // Valgrind did not like the iterator being deleted,
      // thus this complicated loop structure. Did not complain 
      while (it != it_end) 
        {
          if ( (*it)->mCol >= from && (*it)->mCol < to) 
            {
              delete *it;
              PairIterator it2 = it;
              ++it;              
              mPairs.erase(it2);
            }
          else  
            ++it;           
        }

      setChangedLength();

    }

} // namespace alignlib
