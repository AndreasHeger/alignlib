/*
  alignlib - a library for aligning protein sequences

  $Id: ImplFragmentorDiagonals.cpp,v 1.2 2004/01/07 14:35:35 aheger Exp $

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

#include "Fragmentor.h"

#include "Alignment.h"
#include "HelpersAlignment.h"

#include "AlignmentIterator.h"

#include "Alignandum.h"
#include "AlignException.h"

#include "SubstitutionMatrix.h"
#include "HelpersSubstitutionMatrix.h"

#include "Alignator.h"

#include "ImplFragmentorDiagonals.h"

#ifdef WITH_DMALLOC
#include <dmalloc.h>
#endif

using namespace std;

namespace alignlib {

/*---------------------factory functions ---------------------------------- */

  /** make an alignator object, which does a dot-alignment. The default version can be given an AlignmentMatrix-
      object */
Fragmentor * makeFragmentorDiagonals(Score gop, 
				     Score gep, 
				     Alignator * dottor ) {
  
  return new ImplFragmentorDiagonals( gop, gep, gop, gep, dottor );
}


//---------------------------------------------------------------------------------------------------
inline Diagonal calculateDiagonal( const ResiduePAIR & p) { return (p.mCol - p.mRow); }

//----------------------------------------------------------------------------------------------------
  
ImplFragmentorDiagonals::ImplFragmentorDiagonals( Score row_gop, Score row_gep, 
						  Score col_gop, Score col_gep,
						  Alignator * dottor) :
  ImplFragmentor(), 
  mRowGop(row_gop), 
  mRowGep( row_gep), 
  mColGop(col_gop), 
  mColGep( col_gep), 
  mDottor( dottor ) {
}


ImplFragmentorDiagonals::~ImplFragmentorDiagonals() 
{
  debug_func_cerr(5);

}

ImplFragmentorDiagonals::ImplFragmentorDiagonals( const ImplFragmentorDiagonals & src ) : 
   ImplFragmentor(src), 
   mRowGop( src.mRowGop), 
   mRowGep( src.mRowGep), 
   mColGop( src.mColGop), 
   mColGep( src.mColGep), 
   mDottor( src.mDottor)
{
  debug_func_cerr(5);  
}

//------------------------------------------------------------------------------------------
/** CleanUp is empty and does not call base class method, as no shifting is necessary */
void ImplFragmentorDiagonals::cleanUp(const Alignandum * row, const Alignandum *col, Alignment * ali)
{
  debug_func_cerr(5);
}

//------------------------------------------------------------------------------------------
Score ImplFragmentorDiagonals::getGapCost( const ResiduePAIR & p1, const ResiduePAIR & p2 ) const
{

  Position r1 = p1.mRow;
  Position r2 = p2.mRow;  
  Position c1 = p1.mCol;
  Position c2 = p2.mCol;

  Score gap_cost = 0;
  Position d;

  if ((d = (r2 - r1 - 1)) > 0)
      gap_cost += mRowGop + d * mRowGep;

  if ((d = (c2 - c1 - 1)) > 0)
      gap_cost += mColGop + d * mColGep;

  return gap_cost;
}
	
//------------------------------------------------------------------------------------------------
void ImplFragmentorDiagonals::performFragmentation( const Alignandum * row, 
						    const Alignandum * col, 
						    const Alignment * sample) {

    // create dot-matrix (sorted by diagonal)
    Alignment * matrix = makeAlignmentMatrixDiagonal();
    
    // dots are automatically moved to correct position,
    // if only segments are used in row/col
    mDottor->align( row, col, matrix );
    
    // rescore dots (later: use plugin Scoror?), set score to 1
    rescoreAlignment( matrix, 1.0);
    
    Diagonal last_diagonal = MAX_DIAGONAL;
    
    AlignmentIterator it(matrix->begin()), it_end(matrix->end());
    
    Score last_score = 0;
    Alignment * current_ali = makeAlignmentSet();
    
    const ResiduePAIR * last_p = NULL;

    unsigned int length = 0;
    
    for (; it != it_end; ++it) {
	
	const ResiduePAIR * p = &(*it);
	
	Diagonal current_diagonal = calculateDiagonal( *p );
	
	// save alignment if diagonal has been changed
	if (last_diagonal != current_diagonal) {
	    // save fragment with positive score
	    if (last_score > 0 && length > 1) {
		Alignment * new_ali = sample->getNew();
		copyAlignment( new_ali, current_ali);
		new_ali->setScore( last_score );
		mFragments->push_back( new_ali );
	    }
	    
	    // erase alignment
	    current_ali->clear();
	    last_p = NULL;
	    last_diagonal = current_diagonal;
	    last_score = 0;
	    length = 0;
	}
	
	// calculate score along diagonal
	Score new_score;
	if (last_p != NULL)
	    new_score = last_score + p->mScore + getGapCost( *last_p, *p);
	else
	    new_score = last_score + p->mScore;

	// save fragment, if score drops below zero from positive score
	if (new_score <= 0) {
	  // stop existing alignment
	  // save
	  if (last_score > 0 && length > 1) {
	    Alignment * new_ali = sample->getNew();
	    copyAlignment( new_ali, current_ali);
	    new_ali->setScore( last_score );
	    mFragments->push_back( new_ali );
	    length = 0;
	  }

	  current_ali->clear();
	  // start new alignment
	  if (p->mScore > 0) {
	    last_p = p;
	    last_score = p->mScore;
	    current_ali->addPair( new ResiduePAIR( *p ));
	    length = 1;
	  } else {
	    last_p = NULL;
	    last_score = 0;
	    length = 0;
	  }	    
	} else {
	  // continue existing alignment
	  last_score = new_score;
	  last_p = p;
	  current_ali->addPair( new ResiduePAIR( *p ));
	  length++;
	}
    }

    if (last_score > 0 && length > 1) {
	Alignment * new_ali = sample->getNew();
	copyAlignment( new_ali, current_ali);
	new_ali->setScore( last_score );
	mFragments->push_back( new_ali );
    }
    
}
  

} // namespace alignlib



