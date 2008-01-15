/*
  alignlib - a library for aligning protein sequences

  $Id: ImplAlignatorDots.cpp,v 1.2 2004/01/07 14:35:33 aheger Exp $

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

#include "alignlib_fwd.h"
#include "alignlib_interfaces.h"
#include "AlignlibDebug.h"
#include "AlignException.h"
#include "ImplAlignatorDots.h"
#include "Alignandum.h"
#include "ImplAlignmentMatrixRow.h"
#include "Iterator2D.h"
#include "HelpersSubstitutionMatrix.h"

#include "Alignment.h"
#include "HelpersAlignment.h"
#include "HelpersSequence.h"
#include "ImplSequence.h"

using namespace std;

namespace alignlib
{

  /*---------------------factory functions ---------------------------------- */

  //----------------------------------------------------------------------------------------------------------
  /** constructors and destructors */
  ImplAlignatorDots::ImplAlignatorDots( 
		  const HAlignator & dots,
		  Score row_gop, Score row_gep, 
      Score col_gop, Score col_gep ) :
        ImplAlignator(), mDottor (dots),
        mRowGop( row_gop ), mRowGep( row_gep ),
        mColGop( col_gop ), mColGep( col_gep )
        {
          if (mColGop == 0)
            {
              mColGop = mRowGop;
              mColGep = mRowGep;
            }

        }

  //----------------------------------------------------------------------------------------------------------
  ImplAlignatorDots::ImplAlignatorDots( const ImplAlignatorDots & src ) : 
    ImplAlignator( src ), 
    mDottor(src.mDottor) 
    {
      debug_func_cerr(5);

    }

  //----------------------------------------------------------------------------------------------------------
  ImplAlignatorDots::~ImplAlignatorDots()

    {
      debug_func_cerr(5);

    }

  //----------------------------------------------------------------------------------------------------------
  void ImplAlignatorDots::setRowGop( Score gop ) { mRowGop = gop;} 
  void ImplAlignatorDots::setRowGep( Score gep ) { mRowGep = gep;}
  void ImplAlignatorDots::setColGop( Score gop ) { mColGop = gop;}
  void ImplAlignatorDots::setColGep( Score gep ) {mColGep = gep; }

  Score ImplAlignatorDots::getRowGop() { return mRowGop; }
  Score ImplAlignatorDots::getRowGep() { return mRowGep; }
  Score ImplAlignatorDots::getColGop() { return mColGop; }
  Score ImplAlignatorDots::getColGep() { return mColGep; }

  //----------------------------------------------------------------------------------------------------------
  void ImplAlignatorDots::startUp(HAlignment & ali,
		  const HAlignandum & row, 
		  const HAlignandum & col )
    {
      ImplAlignator::startUp(ali, row, col);  
      debug_func_cerr(5);

      mRowLength = mIterator->row_size();
      mColLength = mIterator->col_size();

      mMatrix = makeAlignmentMatrixRow();
      
      // 
      // setup matrix of dots

      // TODO: check with AlignatorPublishDots
      // create dots. This is a patch for AlignatorPublishDots, which
      // does not copy into mMatrix, but returns another pointer.
      mDottor->align( mMatrix, row, col ); 

      // get the number of dots, which corresponds to the length of the
      // alignment in this class. Tell the matrix to sort, etc., at the 
      // same time.
      mNDots = mMatrix->getLength();

      debug_cerr( 5, *mMatrix );

      // get pointers to location of dots(pairs)
      const HAlignment test = makeAlignmentMatrixRow();
      const HImplAlignmentMatrix t = boost::dynamic_pointer_cast< ImplAlignmentMatrix, Alignment>(test); 
      mPairs	= &t->mPairs;
      mRowIndices = t->mIndex;	// these have to be sorted by row, that's why I use AlignmentMatrixRow	

      mTrace   = new int[mNDots];
      mLastDot = -1;
    }

  //-------------------------------------------------------------------------------------------------------
  void ImplAlignatorDots::cleanUp(HAlignment & ali, 
		  const HAlignandum & row, const HAlignandum & col) 
    {
      debug_func_cerr(5);


      if (mTrace != NULL) 
        delete [] mTrace;

      ImplAlignator::cleanUp(ali, row, col );

    }

  //----------------------------------------------------------------------------------------------------------------------------------------
  void ImplAlignatorDots::align(HAlignment & result,
		  const HAlignandum & row, 
		  const HAlignandum & col ) 
    {
      debug_func_cerr(5);

      startUp(result, row, col );

      performAlignment(result, row, col);

      traceBack(result, row, col);

      cleanUp(result, row, col);
    }

  //-----------------------------------------< BackTracke >-------------------------------------------------------------
  void ImplAlignatorDots::traceBack( HAlignment & result,
		  const HAlignandum & row, 
		  const HAlignandum & col ) 
    {
      debug_func_cerr(5);

      int col_res, row_res; 

      int idot   = mLastDot;
      int jleft  = row->getLength();

      while ( idot >= 0) 
        {

          debug_cerr( 5, "-->idot "     << setw(5) << idot      <<
              " col[idot] "  << setw(5) << (*mPairs)[idot]->mCol <<
              " row[idot] "  << setw(5) << (*mPairs)[idot]->mRow <<
              " mTrace[idot] "<< setw(5) << mTrace[idot] );

          row_res = (*mPairs)[idot]->mRow;
          col_res = (*mPairs)[idot]->mCol;

          if (row_res < 0) continue;
          if (col_res < 0) continue;                         
          if (row_res > jleft) break;
          jleft = row_res;                                   // just in case

          result->addPair( new ResiduePAIR(row_res, col_res, (*mPairs)[idot]->mScore) );

          idot = mTrace[idot];
        }

      result->setScore( mScore );
    } 
  //------------------------------------------------------------------------------------------------------------
  // find the index of a residue pairs given row and column
  Position ImplAlignatorDots::getPairIndex( Position r, Position c ) const 
  {
    int x     = mRowIndices[r];     
    bool found = false;  

    if ( x == NO_POS ) 
      return NO_POS;  

    while ((*mPairs)[x]->mRow == r ) {         
      if ((*mPairs)[x]->mCol == c) {             
        found = true;             
        break;         
      }         
      x++;     
    }  

    if (found)         
      return x;     
    else         
      return NO_POS; 
  }

  //----------------------------------------------------------------------------------------------------------------------------------------


  //-----------------------------------------------------------< Alignment subroutine >----------------------------------------------
  void ImplAlignatorDots::performAlignment( HAlignment & ali,
		  const HAlignandum & prow, 
		  const HAlignandum & pcol ) 
    {

      /**
     Overview over the algorithm

     Dots are sorted by row and then by column. This is also the order or processing
     during the algorithm.

     col ->
     row
     |

     ------------|
     |           |
     |           |
     |           |
     |           |
     ------------x

     note: you do not have to consider dots in the same row or column, as it is
     not possible to match the same row to two different columns and vice versa.

     Four different dots are looked up:
     1. The best one in the previous row
     2. The best one in the previous column
     3. The best one directly before the current one
     4. Only if this one is not found, the overall best dot is looked up
     in the upper left part of the matrix( 1 to row-2 and 1 to col-2).
     While proceeding through the same row, the search area only between
     the last lookup and the current col - 1 has to be processed. 

     In order for lookup to be faster, the best dot per column up to then is saved in
     an array bestpercol. This makes the algorithm not strictly Smith-Waterman any more.
     There might be a lower scoring residue, but which is nevertheless closer to the 
     current dot. If two dots in the same col have the same score, the one with the higher
     row is taken. This is ensured, by comparing the score of the dot in the current row
     with the operator >= to the previous best score.

     The array bestpercol has to be updated only after a full row has been processed,
     since its values are needed while the row is still active.



       */

      Score globbest, best;
      Score sa,sb,sc,sd,s;

      int i;

      int last_search_col, xcol = 0;
      Dot best_dot, glob_dot, search_dot, prev_col_dot, prev_row_dot, diagonal_dot, last_dot, current_dot, xdot;

      /* residue numbers of row and column for current dot. These are counted starting at 1 */
      int col_res, row_res;

      int bestpercolstackptr;					/* points to next free element in bestpercolstack */

      /* Allocate Memory for bestpercol and topdot */
      int *bestpercol      = new int[mColLength + 1];
      int *bestpercolstack = new int[mColLength + 1];
      Score *m          = new Score[mNDots];

      best_dot = glob_dot = search_dot = prev_col_dot = prev_row_dot = diagonal_dot = last_dot = -1;
      globbest = best = 0;

      bestpercolstackptr = STACKEMPTY;

      for (col_res = 1; col_res <= mColLength; col_res++) { bestpercol[col_res] = -1; }
      for (i = 0; i < mNDots; i ++) { mTrace[i] = -1; m[i] = 0; }

      last_search_col = 1;
      current_dot = -1;

      // #define DEBUG
      //----------------------------------> main alignment loop <----------------------------------------------------
      for ( current_dot = mRowIndices[1]; current_dot < mNDots; current_dot++ ) {	   /* iterate through nextrow starting at first position */

        if (current_dot < 0) continue;
        row_res = (*mPairs)[current_dot]->mRow;                           /* row_res = row */
        col_res = (*mPairs)[current_dot]->mCol;                           /* col_res = col, wrap around col */

        // some safety checks
#ifdef SAVE
        if (row_res > mRowLength) break;                   /* checking boundaries */
        if (col_res < 1) continue;                /* can be removed, if only dots in boundaries are supplied */
        if (col_res > mColLength) continue;
#endif

#ifdef DEBUG
        std::cout << "--------------------------------------------" << std::endl;
        std::cout << "current_dot=" << current_dot
        << " row_res=" << row_res
        << " col_res=" << col_res
        << " score=" << (*mPairs)[current_dot]->mScore
        << std::endl;
#endif
        /* calculate top row */

        /*------------------------------------------------------------------------------*/
        if ( (last_dot < 0) ||				/* enter first time */
            (row_res  > (*mPairs)[last_dot]->mRow) ) {		/* skip, if not in the same row as last time*/

              /* commit changes to bestpercol from bestpercolstack */
              while( bestpercolstackptr > STACKEMPTY ) {
                xdot = bestpercolstack[--bestpercolstackptr];
                xcol = (*mPairs)[xdot]->mCol;
                if ( (bestpercol[xcol] == -1) || best >= m[bestpercol[xcol]])
                  bestpercol[xcol] = xdot;
              }

#ifdef DEBUG
              for (i = 1; i <= mColLength; i++)
                cout << i << "=" << bestpercol[i] << ";";
              cout << endl;
#endif

              /* init all */
              last_search_col = 1; 
              search_dot = -1; prev_col_dot = -1; prev_row_dot = -1; 
        }

        /*------------------------------------------------------------------------------*/
        /* update prev_row_dot = maximum dot along this row */
        xdot = mRowIndices[row_res-1];
        sc = 0; 
        while ( (xdot > -1)  && 
            ((*mPairs)[xdot]->mRow == row_res-1 ) &&  /* stop, if dot in previous row any more*/
            ((*mPairs)[xdot]->mCol  < col_res-1 )     /* end, if direct contact to new dot*/
        ) { 

          s = m[xdot] + getGapCost( xdot, current_dot);

          if (s > sc) { 
            prev_row_dot = xdot; 
            s = sc; 
          }
          xdot++;
        }

        /*------------------------------------------------------------------------------*/
        /* update prev_col_dot = max scoring dot in previous column */
        if( col_res > 1) {
          prev_col_dot = bestpercol[col_res-1]; 
          sb = m[prev_col_dot] + getGapCost(prev_col_dot, current_dot);
        } else {
          prev_col_dot = -1; 
          sb = 0;
        }

        /*------------------------------------------------------------------------------*/
        /* compute d = match adjacent dot in previous row and column */
        /* look up index in row, col, score for dot; -1 if not found */
        /* diagonal_dot -> dot in previous column, previous row */
        if (col_res > 1) { 
          diagonal_dot = getPairIndex( row_res - 1, col_res - 1);  
          sd = m[diagonal_dot];
        } else {
          diagonal_dot = -1; 
          sd = 0;
        }

        /* update search_dot */
        if ( diagonal_dot < 0 ) { /* only update a if d unoccupied */

          if ( search_dot > -1 )		/* Previous match from last dot in the same current row*/ 
            sa = m[search_dot] + getGapCost( search_dot, current_dot );
          else 
            sa=0; 

          /* Search through area (row_res -2, col_res -2, but start in last_search_col, where
	 we left of the search previously, while being in the same row */

          for (i = last_search_col; i <= col_res-2; i++) {
            xdot = bestpercol[i]; 
            if ( xdot < 0 ) continue;
            s = m[xdot] + getGapCost( xdot, current_dot );
            if (s > sa) { 
              sa = s; 
              search_dot = xdot; 
            }
          }
          last_search_col = col_res - 1;
        } else { 
          search_dot = -1; 
          sa = 0;
        }

        /*------------------------------------------------------------------------------*/
        /*  select best of d|a|b|c */
        best = 0; best_dot = -1; 

        if( sd > best ) { best = sd; best_dot = diagonal_dot; }
        if( sa > best)  { best = sa; best_dot = search_dot; }
        if( sb > best)  { best = sb; best_dot = prev_col_dot; }
        if( sc > best)  { best = sc; best_dot = prev_row_dot; }

        /* record mTraceback */
        best += (*mPairs)[current_dot]->mScore; 

        if (best < 0) { /* local alignment, reset to zero or start new mTrace with single match */
          best    = 0; 
          best_dot = -1;
        }				
        m[current_dot]      = best;
        mTrace[current_dot] = best_dot; 

        if ( best > globbest) { /* save best dot */
          globbest = best; 
          glob_dot  = current_dot;
        }

#ifdef DEBUG
        printf("current_dot %5i; mTrace[current_dot] %5i; m[current_dot] %5.2f; best %5.2f\n",
            current_dot, mTrace[current_dot], m[current_dot], best);
        printf("search_dot %5i; prev_col_dot %5i; prev_row_dot %5i; diagonal_dot %5i\n",
            search_dot,prev_col_dot,prev_row_dot,diagonal_dot);
        printf("sa   %5.2f; sb   %5.2f; sc   %5.2f; sd   %5.2f\n",
            sa, sb, sc, sd);
        printf("last_search_col %5i; glob_dot %5i; globbest %5.2f\n",
            last_search_col, glob_dot, globbest );
#endif

        /* save score as best score per col, if score is higher than previous score */
        if ( (bestpercol[col_res] == -1) || best >= m[bestpercol[col_res]])
          bestpercolstack[bestpercolstackptr++] = current_dot;

        last_dot = current_dot;
      }

      mLastDot= glob_dot;
      mScore  = globbest;

      //--------------> cleaning up <---------------------------------------------------------------

      delete [] bestpercolstack;
      delete [] bestpercol;        
      delete [] m;        

      //--------------------------------------------------------------------------------------------------
    }

} // namespace alignlib
