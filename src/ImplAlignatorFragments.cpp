/*
  alignlib - a library for aligning protein sequences

  $Id: ImplAlignatorFragments.cpp,v 1.3 2004/01/07 14:35:34 aheger Exp $

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
#include <stdio.h>

#include "alignlib.h"
#include "AlignlibDebug.h"
#include "AlignException.h"

#include "Alignandum.h"

#include "ImplAlignatorFragments.h"

#include "Fragmentor.h"

#include "Alignata.h"
#include "HelpersAlignata.h"

#ifdef WITH_DMALLOC
#include <dmalloc.h>
#endif

using namespace std;

namespace alignlib {

#define NOFRAGMENT -1

/*---------------------factory functions ---------------------------------- */

//--------------------------------------------------------------------------------------------------------
  /** constructors and destructors */
  ImplAlignatorFragments::ImplAlignatorFragments( Score row_gop, Score row_gep, 
						  Score col_gop, Score col_gep,
						  Fragmentor * fragmentor) :
    ImplAlignator(), 
    mFragmentor(fragmentor),
    mRowGop( row_gop - row_gep ), mRowGep( row_gop - row_gep ),
    mColGop( col_gop - col_gep ), mColGep( col_gop - col_gep )
  {
    if (mColGop == 0)
      {
	mColGop = mRowGop;
	mColGep = mRowGep;
      }
    
  }
  
  //--------------------------------------------------------------------------------------------------------
  ImplAlignatorFragments::ImplAlignatorFragments( const ImplAlignatorFragments & src ) : 
    ImplAlignator( src ), mFragmentor(src.mFragmentor),
    mRowGop( src.mRowGop), mRowGep( src.mRowGep), mColGop( src.mColGop), mColGep( src.mColGep)
  {
    debug_func_cerr(5);
  }
  
  //--------------------------------------------------------------------------------------------------------
  ImplAlignatorFragments::~ImplAlignatorFragments() 
  {
  debug_func_cerr(5);

  }

  //----------------------------------------------------------------------------------------------------------

  void ImplAlignatorFragments::setRowGop( Score gop ) { mRowGop = gop;} 
  void ImplAlignatorFragments::setRowGep( Score gep ) { mRowGep = gep;}
  void ImplAlignatorFragments::setColGop( Score gop ) { mColGop = gop;}
  void ImplAlignatorFragments::setColGep( Score gep ) {mColGep = gep; }

  Score ImplAlignatorFragments::getRowGop() { return mRowGop; }
  Score ImplAlignatorFragments::getRowGep() { return mRowGep; }
  Score ImplAlignatorFragments::getColGop() { return mColGop; }
  Score ImplAlignatorFragments::getColGep() { return mColGep; }

//--------------------------------------------------------------------------------------------------------
void ImplAlignatorFragments::startUp(const Alignandum * row, const Alignandum *col, Alignata * ali) {
    ImplAlignator::startUp(row, col, ali);  
    debug_func_cerr(5);
    
    mRowLength = mIterator->row_size();
    mColLength = mIterator->col_size();
    
    // setup example alignment type
    Alignata * dummy_alignment = makeAlignataSet();
    
    // create dots
    mFragments = mFragmentor->fragment( row, col, dummy_alignment ); 
    
    // get the number of dots, which corresponds to the length of the
    // alignment in this class. Tell the matrix to sort, etc., at the 
    // same time.
    mNFragments = mFragments->size();

    mTrace        = new int[mNFragments];
    mLastFragment = -1;
}

//--------------------------------------------------------------------------------------------------------
void ImplAlignatorFragments::cleanUp(const Alignandum * row, const Alignandum *col, Alignata * ali) 
{
  debug_func_cerr(5);


  if (mTrace != NULL) 
    delete [] mTrace;

  FragmentVector::iterator it(mFragments->begin()), it_end(mFragments->end());
  for (; it!= it_end; ++it) 
    delete *it;

  delete mFragments;
    
  ImplAlignator::cleanUp(row, col, ali);
  
}

//------------------------------------------------------------------------------------------------------
Alignata * ImplAlignatorFragments::align(const Alignandum * row, const Alignandum * col, Alignata * result) 
{
  debug_func_cerr(5);

  
  startUp(row, col, result);
  
  performAlignment(row, col, result);
    
  traceBack(row, col, result);
  
  cleanUp(row, col, result);
  
  return result;
}

//-----------------------------------------< BackTracke >-----------------------------------------------
void ImplAlignatorFragments::traceBack( const Alignandum * row, const Alignandum * col, Alignata * result) 
{
  debug_func_cerr(5);


    /* trace back along the fragments that are part of
       the optimum trace. Add pairwise alignments to 
       result */

    Fragment ifragment = mLastFragment;

    while ( ifragment >= 0) {
      addAlignata2Alignata( result, (*mFragments)[ifragment] );

#ifdef DEBUG
      cout << "Traceback: fragment " << ifragment << " with score " << (*mFragments)[ifragment]->getScore() << endl;
      writePairAlignment( cout, row, col, (*mFragments)[ifragment] );
#endif

      ifragment = mTrace[ifragment];

    }

    result->setScore( mScore );
} 

//------------------------------------------------------------------------------------------------------



} // namespace alignlib
