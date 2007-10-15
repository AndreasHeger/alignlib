/*
  alignlib - a library for aligning protein sequences

  $Id: ImplAlignataHash.cpp,v 1.4 2004/06/02 12:11:37 aheger Exp $

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
#include "ImplAlignataHash.h"
#include "AlignataIterator.h"

#ifdef WITH_DMALLOC
#include <dmalloc.h>
#endif

using namespace std;

namespace alignlib {

//------------------------------factory functions -----------------------------
Alignata * makeAlignataHash() {
  return new ImplAlignataHash();
}

//------------------------------------< constructors and destructors >-----
ImplAlignataHash::ImplAlignataHash() : ImplAlignata() 
{
  debug_func_cerr(5);

}

ImplAlignataHash::ImplAlignataHash( const ImplAlignataHash& src) : ImplAlignata( src ) 
{
  debug_func_cerr(5);


    mPairs.clear();

    PairIterator it(src.mPairs.begin()), it_end(src.mPairs.end());

    for (;it != it_end; ++it) 
      mPairs.insert( new ResiduePAIR(**it) );
}

ImplAlignataHash::~ImplAlignataHash( ) 
{
  debug_func_cerr(5);

    clear();
}

//------------------------------------------------------------------------------------------------------------
ImplAlignataHash * ImplAlignataHash::getNew() const {
  return new ImplAlignataHash();
}
    
ImplAlignataHash * ImplAlignataHash::getClone() const {
  return new ImplAlignataHash( *this );
}

void ImplAlignataHash::clear() { 
  ImplAlignata::clear();

  PAIRSET::iterator it(mPairs.begin()), it_end(mPairs.end());
  for (;it != it_end; ++it) delete *it;
  mPairs.clear(); 

}

//-----------------------------------------------------------------------------------------------------------   

AlignataConstIterator ImplAlignataHash::begin() const { 
  return AlignataConstIterator( new ImplAlignataHash_ConstIterator(mPairs.begin()));
}

AlignataConstIterator ImplAlignataHash::end()   const { 
  return AlignataConstIterator( new ImplAlignataHash_ConstIterator(mPairs.end())); 
}

AlignataIterator ImplAlignataHash::begin() { 
  return AlignataIterator( new ImplAlignataHash_Iterator(mPairs.begin()));
}

AlignataIterator ImplAlignataHash::end() { 
  return AlignataIterator( new ImplAlignataHash_Iterator(mPairs.end())); 
}


//----------------> accessors <------------------------------------------------------------------------------

Position ImplAlignataHash::getRowFrom() const { return (mPairs.size() > 0) ? (*mPairs.begin())->mRow : 0; }
Position ImplAlignataHash::getColFrom() const { return (mPairs.size() > 0) ? (*mPairs.begin())->mCol : 0; }
Position ImplAlignataHash::getRowTo()   const { return (mPairs.size() > 0) ? (*mPairs.rbegin())->mRow : 0; }
Position ImplAlignataHash::getColTo()   const { return (mPairs.size() > 0) ? (*mPairs.rbegin())->mCol : 0; }
  
ResiduePAIR ImplAlignataHash::getFirstPair() const { return **(mPairs.begin()); }
ResiduePAIR ImplAlignataHash::getLastPair()  const { return **(mPairs.rbegin()); }

void ImplAlignataHash::addPair( ResiduePAIR * new_pair ) 
{
  debug_func_cerr(5);


  if (mPairs.find( new_pair) != mPairs.end()) {
    delete new_pair;
  } else {
    setChangedLength(); 
    mPairs.insert( new_pair ); 
  } 
} 

void ImplAlignataHash::removePair( const ResiduePAIR & old_pair ) 
{
  debug_func_cerr(5);


  setChangedLength(); 
  ResiduePAIR p(old_pair);
  mPairs.erase(&p); 
} 

//----------------------------------------------------------------------------------------------------------
/** retrieves a pair of residues from the alignment */
ResiduePAIR ImplAlignataHash::getPair( Position row ) const {

  ResiduePAIR p(row, 0, 0);
  PairIterator it(mPairs.find( &p ));
  return **it;
} 


} // namespace alignlib
