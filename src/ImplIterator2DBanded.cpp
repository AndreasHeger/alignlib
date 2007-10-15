/**
  alignlib - a library for aligning protein sequences

  $Id: Iterator2D.cpp,v 1.2 2004/01/07 14:35:32 aheger Exp $

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
#include <iomanip>
#include <algorithm>
#include <assert.h>
#include "alignlib.h"
#include "AlignlibDebug.h"
#include "Alignandum.h"
#include "ImplIterator2DBanded.h"

#ifdef WITH_DMALLOC
#include <dmalloc.h>
#endif

using namespace std;

namespace alignlib
{

  /** factory function for creating iterator in band given by offset and width.
      Diagonals are calculated as col - row.

		col_from --------> col_to
      row_from  x  z
      |          x  z
      |         y x  z
      |          y x  z
      |           y x  z  <- mUpperDiagonal
      row_to       <- mLowerDiagonal

      
  */
  Iterator2D * makeIterator2DBanded( const Alignandum * row, const Alignandum * col,
				     const Diagonal lower_diagonal,
				     const Diagonal upper_diagonal )
  {
    return new ImplIterator2DBanded( row, col, lower_diagonal, upper_diagonal );
  }
  
  //--------------------------------------------------------------------------------------
  ImplIterator2DBanded::ImplIterator2DBanded( const Alignandum * row,
					      const Alignandum * col,
					      Diagonal lower_diagonal,
					      Diagonal upper_diagonal ):
    ImplIterator2D( row, col),
    mLowerDiagonal(lower_diagonal), mUpperDiagonal(upper_diagonal)
  {
    debug_func_cerr(5);
    assert(mLowerDiagonal <= mUpperDiagonal);
    
    if (row != NULL && col != NULL)
      resetRanges( row, col);
  }    

  //--------------------------------------------------------------------------------------
  ImplIterator2DBanded::~ImplIterator2DBanded ()
  
{
  debug_func_cerr(5);

  }
  
  //--------------------------------------------------------------------------------------
  ImplIterator2DBanded::ImplIterator2DBanded(const ImplIterator2DBanded & src) :
    ImplIterator2D( src ),
    mLowerDiagonal( src.mLowerDiagonal), mUpperDiagonal( src.mUpperDiagonal )
  
{
  debug_func_cerr(5);

  }

  //--------------------------------------------------------------------------------------
  void ImplIterator2DBanded::resetRanges( const Alignandum * row, const Alignandum * col )
  {
    debug_func_cerr(5);
    
    mRowFrom = std::max( (Position)(row->getFrom()),                  (Position)(col->getFrom() + mUpperDiagonal) );
    mRowTo   = std::min( (Position)(row->getTo()),                    (Position)(col->getTo()   - mUpperDiagonal)  );    
    mColFrom = std::max( (Position)(row->getFrom() + mLowerDiagonal), (Position)(col->getFrom()) );
    mColTo   = std::min( (Position)(row->getTo()   - mLowerDiagonal), (Position)(col->getTo()) );

    debug_cerr( 5, "mRowFrom=" << mRowFrom << " mRowTo=" << mRowTo << " mColFrom=" << mColFrom << " mColTo=" << mColTo );
  }
  
  //--------------------------------------------------------------------------------------  
  /** return a copy of the same iterator
   */
  Iterator2D * ImplIterator2DBanded::getClone() const
    {
      return new ImplIterator2DBanded( *this );
    }
  
  //--------------------------------------------------------------------------------------    
  /** return a new iterator of same type initializes with for row and col
   */
  Iterator2D * ImplIterator2DBanded::getNew( const Alignandum * row, const Alignandum * col ) const
  {
      return new ImplIterator2DBanded( row, col, mLowerDiagonal, mUpperDiagonal );
  }
  
  //--------------------------------------------------------------------------------------
  Iterator2D::const_iterator ImplIterator2DBanded::row_begin ( Position col ) const
  {
    return const_iterator( row_front(col) ) ;
  }
  Iterator2D::const_iterator ImplIterator2DBanded::row_end   ( Position col ) const
  {
    return const_iterator( row_back(col) + 1 ) ;
  }
  
  Iterator2D::const_iterator ImplIterator2DBanded::col_begin ( Position row ) const
  {
    return const_iterator( col_front(row) );
  }
    
  Iterator2D::const_iterator ImplIterator2DBanded::col_end   ( Position row ) const
  {
    return const_iterator( col_back(row) + 1);
  }

  Position ImplIterator2DBanded::row_front ( Position col ) const
  {
    if (col == 0)
      return mRowFrom;
    else
      return ((Position)(col + mLowerDiagonal) >= mRowFrom ) ? (Position)(col + mLowerDiagonal) : 0;
  }

  Position ImplIterator2DBanded::row_back  ( Position col ) const
  {
    if (col == 0)
      return mRowTo;
    else
      return (col + mUpperDiagonal <= mRowTo ) ? (col + mUpperDiagonal) : 0;
  }

  Position ImplIterator2DBanded::col_front ( Position row ) const
  {
    if (row == 0)
      return mColFrom;
    else
      return ( row + mLowerDiagonal >= mColFrom ) ? (row + mLowerDiagonal) : 0;
  }

  Position ImplIterator2DBanded::col_back  ( Position row ) const
  {
    if (row == 0)
      return mColTo;
    else
      return (row + mUpperDiagonal <= mColTo ) ? (row + mUpperDiagonal) : 0;
  }
  
} // namespace alignlib
