/*
  alignlib - a library for aligning protein sequences

  $Id: ImplAlignator.cpp,v 1.3 2005/02/24 11:07:25 aheger Exp $

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
#include "alignlib.h"
#include "AlignlibDebug.h"

#include "Alignator.h"

#include "Alignment.h"
#include "HelpersAlignment.h"

#include "Alignandum.h"
#include "AlignException.h"

#include "Iterator2D.h"
#include "HelpersIterator2D.h"

#include "Scorer.h"
#include "HelpersScorer.h"

#include "ImplAlignator.h"

#ifdef WITH_DMALLOC
#include <dmalloc.h>
#endif

#include <math.h>

using namespace std;

namespace alignlib
{

  //----------------------------------------------------------------------------------------
  ImplAlignator::ImplAlignator()
    {

      mIteratorTemplate = getDefaultIterator2D();
      mIterator = NULL;
      mScorer = NULL;
    }

  ImplAlignator::~ImplAlignator()

    {
      debug_func_cerr(5);

    }

  ImplAlignator::ImplAlignator( const ImplAlignator & src ) : Alignator(src),
  mIterator(src.mIterator),
  mScorer(src.mScorer),
  mIsOwnScorer( src.mIsOwnScorer),
  mIteratorTemplate(src.mIteratorTemplate)

  {
  }

  //-------------------------------------------------------------------------------------------------------------------------------
  void ImplAlignator::setIterator2D( const Iterator2D * iterator ) 
    {
      mIteratorTemplate = iterator;
    }

  void ImplAlignator::setScorer( const Scorer * scorer ) 
    {
      mScorer = scorer;
    }

  void ImplAlignator::startUp( const Alignandum * row, const Alignandum * col, Alignment * ali)

    {
      debug_func_cerr(5);

      row->prepare();
      col->prepare();

      debug_cerr( 5, "starting alignment for row=" << row->getFrom() << "-" << row->getTo() 
          << " col=" << col->getFrom() << "-" << col->getTo() );
      
      mRowLength = row->getLength();

      mIterator = mIteratorTemplate->getNew( row, col );
      
      debug_cerr( 5, "setting iterator to ranges: row=" 
          << *mIterator->row_begin() << "-" <<  *mIterator->row_end() << ":" << mIterator->row_size() << " col=" 
          << *mIterator->col_begin() << "-" <<  *mIterator->col_end() << ":" << mIterator->col_size() );
      
      
      if (mScorer == NULL)
        {
          mScorer = makeScorer( row, col);
          mIsOwnScorer = true;
        }
      else
        {
          mScorer = mScorer->getNew( row, col );
          mIsOwnScorer = true;
        }

      ali->clear();
    }

  void ImplAlignator::cleanUp(const Alignandum * row, const Alignandum *col, Alignment * ali) 
    {
      debug_func_cerr(5);


      if (mIsOwnScorer)
        {
          delete mScorer;
          mScorer = NULL;
          mIsOwnScorer = false;
        }

      delete mIterator;
      mIterator = NULL;

      /* round score to integer. This will get rid
       of ???
       */
      ali->setScore( round(ali->getScore()) );

    }

} // namespace alignlib
