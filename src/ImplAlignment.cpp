/*
  alignlib - a library for aligning protein sequences

  $Id: ImplAlignment.cpp,v 1.4 2004/06/02 12:11:37 aheger Exp $

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
#include "ImplAlignment.h"
#include "Alignandum.h"
#include "Alignatum.h"
#include "alignlib.h"
#include "AlignlibDebug.h"
#include "AlignmentIterator.h"
#include "Alignment.h"
#include "AlignException.h"

using namespace std;

namespace alignlib 
{

  //--------------------------------------------------------------------------------------------------------------
  // constructors and desctructors
  //--------------------------------------------------------------------------------------------------------------
  ImplAlignment::ImplAlignment() : mChangedLength(true), mLength(0), mScore(0), mNumGaps( 0 ), 
  mRowFrom(NO_POS),
  mRowTo(NO_POS),
  mColFrom( NO_POS),
  mColTo(NO_POS)
  {
  }

  ImplAlignment::ImplAlignment( const ImplAlignment & src ) :
    mChangedLength(src.mChangedLength), 
    mLength( src.mLength), 
    mScore (src.mScore), 
    mNumGaps (src.mNumGaps),
    mRowFrom( src.mRowFrom ),
    mRowTo( src.mRowTo ),
    mColFrom( src.mColFrom ),
    mColTo( src.mColTo )
    {
      debug_func_cerr(5);

    }

  ImplAlignment::~ImplAlignment() 
    {
      debug_func_cerr(5);

    }

  //--------------------------------------------------------------------------------------------------------------
  void ImplAlignment::setChangedLength() 
  {
    mChangedLength = true;
  }
  //--------------------------------------------------------------------------------------------------------------
  void ImplAlignment::setScore( Score score ) 
  {
    mScore = score;
  }
  
  //--------------------------------------------------------------------------------------------------------------  
  Position ImplAlignment::getRowFrom() const { return mRowFrom; }
  Position ImplAlignment::getRowTo() const { return mRowTo; }
  Position ImplAlignment::getColFrom() const { return mColFrom; }
  Position ImplAlignment::getColTo() const { return mColTo; }
  
  //--------------------------------------------------------------------------------------------------------------
  Score ImplAlignment::getScore() const 
  {
    return mScore;
  }
  
  //--------------------------------------------------------------------------------------------------------------
  Position ImplAlignment::getLength() const 
  {

    if (mChangedLength)
      calculateLength();

    return mLength;
  }

  //--------------------------------------------------------------------------------------------------------------
  void ImplAlignment::clear() 
    {
      debug_func_cerr(5);

      mChangedLength = false;
      mLength = 0;
      mScore = 0;
      mNumGaps = 0;
    }

  Position ImplAlignment::getNumAligned() const
  {	
	  if (mChangedLength)
		  calculateLength();
	  return mLength - mNumGaps;
  }
  
  //--------------------------------------------------------------------------------------------------------------
  Position ImplAlignment::getNumGaps() const 
  {

    return mNumGaps;
  }
  //--------------------------------------------------------------------------------------------------------------
  bool ImplAlignment::isEmpty() const 
  {
    return (getLength() == 0);
  }

  //-----------------------------------------------------------------------------------------------------------   
  void ImplAlignment::write( std::ostream& output ) const 
  {
    debug_func_cerr(5);


    output << "Length: " << getLength() << "\tScore: " << getScore() << "\tGaps: " << getNumGaps() << endl;
    output << "Row\tColumn\tScore\t" << endl;

    AlignmentConstIterator it(begin());
    AlignmentConstIterator it_end(end());

    for (; it != it_end; ++it)
      output << *it << endl;

  }
  
  //--------------------------------------------------------------------------------------------------------------
  /** Read data members from stream */
  void ImplAlignment::read( std::istream & input) 
  {
  }

  //--------------------------------------------------------------------------------------------------------------
  void ImplAlignment::moveAlignment( Position row_offset, Position col_offset) 
    {
      debug_func_cerr(5);

      if (mLength == 0) return;
      
      AlignmentIterator it     = begin();
      AlignmentIterator it_end = end();

      for (;it != it_end; ++it) 
      {
        ResiduePAIR & p = (*it);
        p.mRow += row_offset;
        p.mCol += col_offset;
      }
      
      mRowFrom += row_offset;
      mRowTo += row_offset;
      mColFrom += col_offset;
      mColTo += col_offset;
      
    }

  //--------------------------------------------------------------------------------------------------------------
  void ImplAlignment::setLength( Position length ) const 
  {
    mLength = length;
  }

  //--------------------------------------------------------------------------------------------------------------
  void ImplAlignment::setNumGaps( Position num_gaps ) const 
  {
    mNumGaps = num_gaps;
  }

  //--------------------------------------------------------------------------------------------------------------  
  void ImplAlignment::addPair( ResiduePAIR * pair )
  {
	  Position row = pair->mRow;
	  Position col = pair->mCol;
	  
	  if (mRowFrom == NO_POS)
	  {
		  mRowFrom = row;
		  mColFrom = col;
		  mRowTo = row + 1;
		  mColTo = col + 1;
	  }
	  else
	  {
		  if (row < mRowFrom) 
			  mRowFrom = row;
		  else if ( ++row > mRowTo )
			  mRowTo = row;
		  
		  if (col < mColFrom) 
			  mColFrom = col;
		  else if (++col > mColTo)   
			  mColTo = col;
	  }
		  
  }

  //--------------------------------------------------------------------------------------------------------------  
  void ImplAlignment::removePair( const ResiduePAIR & pair )
  {
	  debug_func_cerr( 5 );
	  if (pair.mRow == mRowFrom || pair.mRow == mRowTo || pair.mCol == mColFrom || pair.mCol == mColTo )
		  updateBoundaries();
	  setChangedLength();
  }
  
  //--------------------------------------------------------------------------------------------------------------
  void ImplAlignment::addPair( Position row, Position col, Score score ) 
    {
      debug_func_cerr(5);

      addPair(new ResiduePAIR( row,col,score) );
    } 

  //--------------------------------------------------------------------------------------------------------------
  void ImplAlignment::calculateLength() const 
  {
    debug_func_cerr(5);

    AlignmentConstIterator it(begin());
    AlignmentConstIterator it_end(end());

    mLength = 0;
    mNumGaps = 0;
    
    if (it == it_end)
    {
    	mRowFrom = mRowTo = mColFrom = mColTo = NO_POS;
    	return;
    }
    
    mRowFrom = it->mRow;
    mColFrom = it->mCol;
    mRowTo = it->mRow;
    mColTo = it->mCol;

    ++it;
    Position row_last = mRowFrom;
    Position col_last = mColFrom;
    Position d;
    ++mLength;
    
    for (; it != it_end; ++it) 
    {
    	Position row = it->mRow;
    	Position col = it->mCol;
    	
		// get maximum boundaries
    	if (row < mRowFrom) mRowFrom = row;
    	if (col < mColFrom) mColFrom = col;
    	if (row > mRowTo)   mRowTo = row;
    	if (col > mColTo)   mColTo = col;
	
    	++mLength;
    	if ( (d = row - row_last - 1) > 0 ) 
    	{
    		mLength += d; mNumGaps += d;
    	}

    	if ( (d = col - col_last - 1) > 0 ) 
    	{
    		mLength += d; mNumGaps += d;
    	}

    	row_last = row;
    	col_last = col;
    }

    ++mRowTo;
    ++mColTo;
    mChangedLength = false;
  }  

  //-----------------------------------------------------------------------------------------------------------   
  /** switch row and column in the alignment. Use more efficient implementations in derived classes. */
  void ImplAlignment::switchRowCol() 
    {

      debug_func_cerr(5);

      Alignment * copy = getClone();

      AlignmentIterator it     = copy->begin();
      AlignmentIterator it_end = copy->end();

      clear();
      for (;it != it_end; ++it) 
        {
          // copy over residue pairs from copy reversing row and column
          addPair( new ResiduePAIR( it->mCol, it->mRow, it->mScore ) );
        }	

      delete copy;

      std::swap( mRowFrom, mColFrom );
      std::swap( mRowTo, mColTo );
      
      return;
    }

  //-----------------------------------------------------------------------------------------------------------   
  //--------------> mapping functions <----------------------------------------------------------------------------
  /** default implementation for mapping
   * 
   * This function iterates through the alignment and is
   * very time inefficient and ignores the search option.
   */
  Position ImplAlignment::mapRowToCol( Position pos, SearchType search ) const 
  {
    debug_func_cerr(5);

    if (getLength() == 0)
      return NO_POS;

    AlignmentConstIterator it = begin();
    AlignmentConstIterator it_end = end();

    for (; it != it_end; ++it) 
      {
        if ((*it).mRow == pos)
          return (*it).mCol;
      }     

    return NO_POS;
  }

  //-----------------------------------------------------------------------------------------------------------   
  /** default implementation for mapping
   * 
   * This function iterates through the alignment and is
   * very time inefficient and ignores the search option.
   */
  Position ImplAlignment::mapColToRow( Position pos, SearchType search ) const 
  {

    if (getLength() == 0)
      return NO_POS;

    AlignmentConstIterator it = begin();
    AlignmentConstIterator it_end = end();

    while (it != it_end) {
      if ((*it).mCol == pos)
        return (*it).mRow;
      ++it;
    }

    return NO_POS;
  }

  //-----------------------------------------------------------------------------------------------------------   
  /** This is a generic routine. It creates a new alignment by making a copy of the old one.
   */
  void ImplAlignment::removeRowRegion( Position from, Position to) 
  {

    const Alignment * copy = getClone();  

    AlignmentConstIterator it     = copy->begin();
    AlignmentConstIterator it_end = copy->end();

    clear();

    mScore = copy->getScore();

    for (; it != it_end; ++it) 
      {
        if ( (*it).mRow < from || (*it).mRow >= to)
          addPair( new ResiduePAIR(*it) );
      }

    delete copy;

    updateBoundaries();
    setChangedLength();
    return;
  }

  //-----------------------------------------------------------------------------------------------------------   
  /** This is a generic routine. It creates a new alignment by making a copy of the old one.
   */
  void ImplAlignment::removeColRegion( Position from, Position to) 
  {

    const Alignment * copy = getClone();  

    AlignmentConstIterator it     = copy->begin();
    AlignmentConstIterator it_end = copy->end();

    clear();

    mScore = copy->getScore();

    for (; it != it_end; ++it) 
      {
        const ResiduePAIR & p = (*it);
        if (p.mCol < from || p.mCol >= to)
          addPair( new ResiduePAIR(p) );
      }     

    delete copy;

    updateBoundaries();
    setChangedLength();
    return;
  }


  //--------------------------------------------------------------------------------------------------------------



} // namespace alignlib
