/*
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

#include "alignlib.h"
#include "AlignlibDebug.h"
#include "Alignandum.h"
#include "ImplIterator2DFull.h"

#ifdef WITH_DMALLOC
#include <dmalloc.h>
#endif

using namespace std;

namespace alignlib
{

  // factory function for creating iterator over full matrix
  Iterator2D * makeIterator2DFull( const Alignandum * row, const Alignandum * col )
  {
    return new ImplIterator2DFull( row, col );
  }
  
  //--------------------------------------------------------------------------------------
  ImplIterator2DFull::ImplIterator2DFull( const Alignandum * row, const Alignandum * col) :
    ImplIterator2D( row, col )
  {
    debug_func_cerr(5);
  }    
  
  //--------------------------------------------------------------------------------------
  ImplIterator2DFull::~ImplIterator2DFull ()
  {
    debug_func_cerr(5);
  }
  
  //--------------------------------------------------------------------------------------
  ImplIterator2DFull::ImplIterator2DFull(const ImplIterator2DFull & src) :
    ImplIterator2D( src )
  
{
  debug_func_cerr(5);

  }

  //--------------------------------------------------------------------------------------  
  /** return a copy of the same iterator
   */
  Iterator2D * ImplIterator2DFull::getClone() const
    {
      return new ImplIterator2DFull( *this );
    }

  //--------------------------------------------------------------------------------------
  void ImplIterator2DFull::resetRanges( const Alignandum * row, const Alignandum * col )
  {
    ImplIterator2D::resetRanges( row, col);
  }
  
  
  //--------------------------------------------------------------------------------------    
  /** return a new iterator of same type initializes with for row and col
   */
  Iterator2D * ImplIterator2DFull::getNew( const Alignandum * row, const Alignandum * col ) const
    {
      return new ImplIterator2DFull( row, col );
    }
  
  //--------------------------------------------------------------------------------------
  Iterator2D::const_iterator ImplIterator2DFull::row_begin ( Position col ) const
  {
    return const_iterator( mRowFrom );
  }
  Iterator2D::const_iterator ImplIterator2DFull::row_end   ( Position col ) const
  {
    return const_iterator( mRowTo );
  }
  
  Iterator2D::const_iterator ImplIterator2DFull::col_begin ( Position row ) const
  {
    return const_iterator( mColFrom );
  }
    
  Iterator2D::const_iterator ImplIterator2DFull::col_end   ( Position row ) const
  {
    return const_iterator( mColTo );
  }

  Position ImplIterator2DFull::row_front ( Position col ) const
  {
    return mRowFrom;
  }

  Position ImplIterator2DFull::row_back  ( Position col ) const
  {
    return mRowTo;
  }

  Position ImplIterator2DFull::col_front ( Position row ) const
  {
    return mColFrom;
  }

  Position ImplIterator2DFull::col_back  ( Position row ) const
  {
    return mColTo;
  }
  Position ImplIterator2DFull::row_size ( Position col ) const
  {
    return mRowTo - mRowFrom;
  }

  Position ImplIterator2DFull::col_size ( Position row ) const
  {
    return mColTo - mColFrom;
  }
  
} // namespace alignlib
