/*
  alignlib - a library for aligning protein sequences

  $Id: ImplAlignatorDiagonalDotsSquared.cpp,v 1.2 2004/01/07 14:35:33 aheger Exp $

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

#include <map>
#include <vector>

#include "alignlib.h"
#include "AlignlibDebug.h"
#include "AlignException.h"
#include "ImplAlignatorDiagonalDotsSquared.h"
#include "Alignandum.h"
#include "ImplAlignmentMatrixRow.h"

#include "HelpersSubstitutionMatrix.h"

#include "Alignment.h"
#include "HelpersAlignment.h"


#ifdef WITH_DMALLOC
#include <dmalloc.h>
#endif

using namespace std;

namespace alignlib {

#define NODOT -1

/*---------------------factory functions ---------------------------------- */

    /** make an alignator object, which does a dot-alignment. The default version can be given an AlignmentMatrix-
	object */
Alignator * makeAlignatorDiagonalDotsSquared(Score gop, Score gep, Alignator * alignator, const SubstitutionMatrix * subst_matrix) {

  if (!subst_matrix) 
      return new ImplAlignatorDiagonalDotsSquared( getDefaultSubstitutionMatrix(), 
				    gop, gep, gop, gep, alignator );
  else
      return new ImplAlignatorDiagonalDotsSquared( subst_matrix, 
				    gop, gep, gop, gep, alignator );
  
}

//----------------------------------------------------------------------------------------------------------------------------------------
    /** constructors and destructors */
ImplAlignatorDiagonalDotsSquared::ImplAlignatorDiagonalDotsSquared( const SubstitutionMatrix * subst_matrix,
				      Score row_gop, Score row_gep, 
				      Score col_gop, Score col_gep,
				      Alignator * dots) :
    ImplAlignatorDots( subst_matrix, row_gop - row_gep, row_gep, col_gop - col_gep, col_gep, dots) {
}

//------------------------------------------------------------------------------------------
ImplAlignatorDiagonalDotsSquared::ImplAlignatorDiagonalDotsSquared( const ImplAlignatorDiagonalDotsSquared & src ) : 
  ImplAlignatorDots( src ) 
{
  debug_func_cerr(5);

}

//------------------------------------------------------------------------------------------
ImplAlignatorDiagonalDotsSquared::~ImplAlignatorDiagonalDotsSquared() 
{
  debug_func_cerr(5);

}

ImplAlignatorDiagonalDotsSquared * ImplAlignatorDiagonalDotsSquared::getClone() const 
{
  return new ImplAlignatorDiagonalDotsSquared( *this );
}


//------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------
Score ImplAlignatorDiagonalDotsSquared::getGapCost( Dot x1, Dot x2 ) const {

  Position c1 = (*mPairs)[x1]->mCol;
  Position c2 = (*mPairs)[x2]->mCol;  
  Position r1 = (*mPairs)[x1]->mRow;
  Position r2 = (*mPairs)[x2]->mRow;  

  Score gap_cost = 0;
  Position d;

  if ((d = (r2 - r1)) > 1)
      gap_cost += mRowGop + d * mRowGep;

  if ((d = (c2 - c1)) > 1)
      gap_cost += mColGop + d * mColGep;

  return gap_cost;
}

//-------------------------------------------< Alignment subroutine >----------------------------------------------
void ImplAlignatorDiagonalDotsSquared::performAlignment( const Alignandum * prow, const Alignandum * pcol, Alignment * ali) {

  /**
     Overview over the algorithm
     
     1. Dots are sorted by row and then by column

     col ->
     row
     |
     
     ------------|
     |           |
     |           |
     |           |
 c-> |           |
     ------------x

     note: you do not have to consider dots in the same row or column, as it is
     not possible to match the same row to two different columns and vice versa.
     
     

  */
//  #define DEBUG
//  #define DEBUG2

  Dot global_best_dot = NODOT;
  Score global_best_score = 0;

  // create data structure for search region
  // sort dots in search region by increasing column
  typedef multimap <Diagonal, Dot> MyDotSet;

  MyDotSet search_region;
  
  // array with scores of dots
  vector<Score> scores(mNDots,0);

  // array with dots in current row
  vector<Dot> dot_stack(mRowLength, NODOT);

  unsigned int num_row_dots = 0;
  Position last_diagonal = MAX_DIAGONAL;

  //----------------------------------> main alignment loop <----------------------------------------------------
  for ( Dot current_dot = 0; current_dot < mNDots; current_dot++ ) {	   

    Position current_row = (*mPairs)[current_dot]->mRow;                         
    Position current_col = (*mPairs)[current_dot]->mCol;                         

    /* if a new row is entered, enter dots from stack to search-area */
    if (current_row != last_row) {
      while (num_row_dots > 0) {
	Dot dot = dot_stack[--num_row_dots];
	search_region.insert(pair<Position, Dot>((*mPairs)[dot]->mCol, dot));
      }
      last_row = current_row;
    }

    /* search search-area: always lookup starting at col 1 until current_col - 1. Try to find
       a positive trace leading to current dot. If it were negative, it would not be part of
       the optimum alignment up to current_dot. */
    Dot search_best_dot   = NODOT;
    Score search_best_score = 0;


#ifdef DEBUG2
    cout << "SEARCH_AREA" << endl;
    for (MyDotSet::iterator it = search_region.begin(); it != search_region.end(); ++it)
      cout << "  [" << (*it).first << ", " << (*it).second << "]" << endl;
#endif
    
    MyDotSet::const_iterator it(search_region.begin()), it_end(search_region.end());
    while (it != it_end && ((*it).first < current_col)) {
      
      Dot const search_dot = (*it).second;
      Score search_score = scores[search_dot];
      
      if (search_score > 0) {
	search_score += getGapCost( search_dot, current_dot);
	if (search_score >= search_best_score) {
	  search_best_score = search_score;
	  search_best_dot   = search_dot;
	}
      }

      ++it;
    }

    /* no positive trace found, new trace starts at current dot */
    if (search_best_dot == NODOT)
      search_best_score = (*mPairs)[current_dot]->mScore;
    else
      search_best_score += (*mPairs)[current_dot]->mScore;

#ifdef DEBUG
    cout << "current_dot=" << current_dot << " current_row=" << current_row << " current_col=" << current_col << endl;
    cout << "search_best_dot=" << search_best_dot << " search_best_score=" << search_best_score << endl;
#endif
    
    /* do local alignment, traces with score <= 0 are skipped */
    if (search_best_score < 0)
      continue;
    
    scores[current_dot] = search_best_score;
    mTrace[current_dot] = search_best_dot;

    dot_stack[num_row_dots++] = current_dot;    

    /* remember end point of best trace */
    if (search_best_score > global_best_score) {
      global_best_score = search_best_score;
      global_best_dot   = current_dot;
    }
    
  } /* end of alignment loop */

  mLastDot= global_best_dot;
  mScore  = global_best_score;

#ifdef DEBUG
  cout << "global_best_dot=" << global_best_dot << " global_best_score=" << global_best_score << endl;
#endif

  //--------------> cleaning up <---------------------------------------------------------------
  //--------------------------------------------------------------------------------------------------
}

} // namespace alignlib







