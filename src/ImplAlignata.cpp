/*
  alignlib - a library for aligning protein sequences

  $Id: ImplAlignata.cpp,v 1.4 2004/06/02 12:11:37 aheger Exp $

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
#include <string>
#include <stdio.h>
#include "ImplAlignata.h"
#include "Alignandum.h"
#include "Alignatum.h"
#include "alignlib.h"
#include "AlignlibDebug.h"
#include "AlignataIterator.h"
#include "Alignata.h"
#include "AlignException.h"

#ifdef WITH_DMALLOC
#include <dmalloc.h>
#endif

using namespace std;

namespace alignlib {

//--------------------------------------------------------------------------------------------------------------
// constructors and desctructors
//--------------------------------------------------------------------------------------------------------------
ImplAlignata::ImplAlignata() : mChangedLength(true), mLength(0), mScore(0), mNumGaps( 0 ) {
}

ImplAlignata::ImplAlignata( const ImplAlignata & src ) :
  mChangedLength(src.mChangedLength), 
  mLength( src.mLength), 
  mScore (src.mScore), 
  mNumGaps (src.mNumGaps) 
{
  debug_func_cerr(5);

}

ImplAlignata::~ImplAlignata() 
{
  debug_func_cerr(5);

}

//--------------------------------------------------------------------------------------------------------------
void ImplAlignata::setChangedLength() {
  mChangedLength = true;
}
//--------------------------------------------------------------------------------------------------------------
void ImplAlignata::setScore( Score score ) {
    mScore = score;
}
//--------------------------------------------------------------------------------------------------------------
Score ImplAlignata::getScore() const {
    return mScore;
}
//--------------------------------------------------------------------------------------------------------------
Position ImplAlignata::getLength() const {

  if (mChangedLength)
    calculateLength();

  return mLength;
}
//--------------------------------------------------------------------------------------------------------------
void ImplAlignata::clear() 
{
  debug_func_cerr(5);

    mChangedLength = false;
    mLength = 0;
    mScore = 0;
    mNumGaps = 0;
}

//--------------------------------------------------------------------------------------------------------------
Position ImplAlignata::getNumGaps() const {
  if (mChangedLength)
    calculateLength();

  return mNumGaps;
}
//--------------------------------------------------------------------------------------------------------------
bool ImplAlignata::isEmpty() const {
    return (getLength() == 0);
}

//-----------------------------------------------------------------------------------------------------------   
void ImplAlignata::write( std::ostream& output ) const 
{
  debug_func_cerr(5);


    output << "Length: " << getLength() << "\tScore: " << getScore() << "\tGaps: " << getNumGaps() << endl;
    output << "Row\tColumn\tScore\t" << endl;

    AlignataConstIterator it(begin());
    AlignataConstIterator it_end(end());
    
    for (; it != it_end; ++it)
      output << *it << endl;

}

//--------------------------------------------------------------------------------------------------------------
/** Read data members from stream */
void ImplAlignata::read( std::istream & input) {
}

//--------------------------------------------------------------------------------------------------------------
void ImplAlignata::moveAlignment( Position row_offset, Position col_offset) 
{
  debug_func_cerr(5);

    
  AlignataIterator it     = begin();
  AlignataIterator it_end = end();
  
  for (;it != it_end; ++it) {
    ResiduePAIR & p = (*it);
    p.mRow += row_offset;
    p.mCol += col_offset;
  }
  
}
   
//--------------------------------------------------------------------------------------------------------------
void ImplAlignata::setLength( Position length ) const {
    mLength = length;
}

//--------------------------------------------------------------------------------------------------------------
void ImplAlignata::setNumGaps( Position num_gaps ) const {
    mNumGaps = num_gaps;
}

//--------------------------------------------------------------------------------------------------------------
void ImplAlignata::addPair( Position row, Position col, Score score ) 
{
  debug_func_cerr(5);

  addPair(new ResiduePAIR( row,col,score) );
} 

//--------------------------------------------------------------------------------------------------------------
void ImplAlignata::calculateLength() const 
{
  debug_func_cerr(5);


  Position row_last = getRowFrom();
  Position col_last = getColFrom();
  Position row_cur, col_cur;

  mLength = 0;
  mNumGaps = 0;
  mChangedLength = false;

  // this is a patch, there is something wrong with empty alignments,
  // for vectors (and sets)
  
  if (row_last == NO_POS || col_last == NO_POS)
      return;

  AlignataConstIterator it(begin());
  AlignataConstIterator it_end(end());

  Position d;
  for (; it != it_end; ++it) {
      ResiduePAIR p(*it);
  
      row_cur = p.mRow;
      col_cur = p.mCol;
    
      mLength++;
      if ( (d = row_cur - row_last - 1) > 0 ) {
	  mLength += d; mNumGaps += d;
      }

      if ( (d = col_cur - col_last - 1) > 0 ) {
	  mLength += d; mNumGaps += d;
      }
    
      row_last = row_cur;
      col_last = col_cur;
  }

}  

//-----------------------------------------------------------------------------------------------------------   
/** switch row and column in the alignment. Use more efficient implementations in derived classes. */
void ImplAlignata::switchRowCol() {

  Alignata * copy = getClone();
  
  AlignataIterator it     = copy->begin();
  AlignataIterator it_end = copy->end();
  
  clear();
  for (;it != it_end; ++it) 
  	{
	  // copy over residue pairs from copy reversing row and column
	  addPair( new ResiduePAIR( it->mCol, it->mRow, it->mScore ) );
  	}	
  
  delete copy;

  setChangedLength();
  return;
}
  
//-----------------------------------------------------------------------------------------------------------   
//--------------> mapping functions <----------------------------------------------------------------------------
/** default implementation for mapping
 * 
 * This function iterates through the alignment and is
 * very time inefficient and ignores the search option.
 */
Position ImplAlignata::mapRowToCol( Position pos, SearchType search ) const {
  
  if (getLength() == 0)
    return NO_POS;

  AlignataConstIterator it = begin();
  AlignataConstIterator it_end = end();
  
  while (it != it_end) {
    if ((*it).mRow == pos)
      return (*it).mCol;
    ++it;
  }
  
  return NO_POS;
}

//-----------------------------------------------------------------------------------------------------------   
/** default implementation for mapping
 * 
 * This function iterates through the alignment and is
 * very time inefficient and ignores the search option.
 */
Position ImplAlignata::mapColToRow( Position pos, SearchType search ) const {

  if (getLength() == 0)
    return NO_POS;

  AlignataConstIterator it = begin();
  AlignataConstIterator it_end = end();
  
  while (it != it_end) {
    if ((*it).mCol == pos)
      return (*it).mRow;
    ++it;
  }
  
  return 0;
}

//-----------------------------------------------------------------------------------------------------------   
/** This is a generic routine. It creates a new alignment by making a copy of the old one.
*/
void ImplAlignata::removeRowRegion( Position from, Position to) {

  const Alignata * copy = getClone();  
  
  AlignataConstIterator it     = copy->begin();
  AlignataConstIterator it_end = copy->end();

  clear();

  mScore = copy->getScore();

  for (; it != it_end; ++it) {
    if ( (*it).mRow < from || (*it).mRow >= to)
      addPair( new ResiduePAIR(*it) );
  }
  
  delete copy;

  setChangedLength();
  return;
}

//-----------------------------------------------------------------------------------------------------------   
/** This is a generic routine. It creates a new alignment by making a copy of the old one.
*/
void ImplAlignata::removeColRegion( Position from, Position to) {

  const Alignata * copy = getClone();  
  
  AlignataConstIterator it     = copy->begin();
  AlignataConstIterator it_end = copy->end();

  clear();

  mScore = copy->getScore();
  
  for (; it != it_end; ++it) {
    const ResiduePAIR & p = (*it);
    if (p.mCol < from || p.mCol >= to)
      addPair( new ResiduePAIR(p) );
  }
  
  delete copy;

  setChangedLength();
  return;
}


// //--------------------------------------------------------------------------------------------------------------



} // namespace alignlib
