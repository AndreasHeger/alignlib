/*
  alignlib - a library for aligning protein sequences

  $Id: ImplAlignatorDotsNew.cpp,v 1.2 2004/01/07 14:35:34 aheger Exp $

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
#include "AlignlibException.h"
#include "ImplAlignatorDots.h"
#include "Alignandum.h"
#include "ImplAlignmentMatrixRow.h"

#include "HelpersSubstitutionMatrix.h"

#include "Alignment.h"
#include "HelpersAlignment.h"

#ifdef DEBUG
#include "stdio.h"
#endif


#ifdef WITH_DMALLOC
#include <dmalloc.h>
#endif

using namespace std;

namespace alignlib {

#define NODOT -1

  /*---------------------factory functions ---------------------------------- */

  /** make an alignator object, which does a dot-alignment. The default version can be given an AlignmentMatrix-
	object */
  Alignator * makeAlignatorDots(Score gop, Score gep, Alignator * alignator, const SubstitutionMatrix * subst_matrix) 
  {
    if (!subst_matrix) 
      return new ImplAlignatorDots( getDefaultSubstitutionMatrix(), 
          gop, gep, gop, gep, alignator );
    else
      return new ImplAlignatorDots( subst_matrix, 
          gop, gep, gop, gep, alignator );

  }

  //----------------------------------------------------------------------------------------------------------------------------------------
  /** constructors and destructors */
  ImplAlignatorDots::ImplAlignatorDots( const SubstitutionMatrix * subst_matrix,
      Score row_gop, Score row_gep, 
      Score col_gop, Score col_gep,
      Alignator * dots) :
        ImplAlignator( subst_matrix, row_gop - row_gep, row_gep, col_gop - col_gep, col_gep), mDottor (dots) 
        {
        }

  //----------------------------------------------------------------------------------------------------------------------------------------
  ImplAlignatorDots::ImplAlignatorDots( const ImplAlignatorDots & src ) : ImplAlignator( src ), mDottor(src.mDottor) 
  {
    debug_func_cerr(5);

  }

  //----------------------------------------------------------------------------------------------------------------------------------------
  ImplAlignatorDots::~ImplAlignatorDots() 
    {
      debug_func_cerr(5);

    }

   //----------------------------------------------------------------------------------------------------------
  ImplAlignatorDotsNew * ImplAlignatorDotsNew::getClone() const 
  {
    return new ImplAlignatorDotsNew( *this );
  }

  //----------------------------------------------------------------------------------------------------------------------------------------
  void ImplAlignatorDots::startUp(const HAlignandum row, const HAlignandumcol, HAlignment ali) 
    {
    ImplAlignator::startUp(row, col, ali);  

    debug_cerr( 5, "void AlignatorDots::Initialize()" );

    // setup matrix of dots
    mMatrix = (ImplAlignmentMatrixRow*)makeAlignmentMatrixRow();

    // create dots
    mDottor->align( row, col, mMatrix ); 

    // get the number of dots, which corresponds to the length of the
    // alignment in this class. Tell the matrix to sort, etc., at the 
    // same time.
    mNDots = mMatrix->getLength();

    debug_cerr( 5, "matrix=" << *mMatrix );

    // setup pointers to location of dots(pairs)
    mPairs	= mMatrix->mPairs;
    mRowIndices = mMatrix->mIndex;	// these have to be sorted by row, that's why I use AlignmentMatrixRow	

    mTrace   = new int[mNDots];
    mLastDot = -1;
  }

  //-----------------------------------------------------------------------------------------------------------------------------
  void ImplAlignatorDots::cleanUp(const HAlignandum row, const HAlignandumcol, HAlignment ali) 
    {
      debug_func_cerr(5);


      if (mTrace != NULL) 
        delete [] mTrace;

      delete mMatrix;

      ImplAlignator::cleanUp(row, col, ali);

    }

  //----------------------------------------------------------------------------------------------------------------------------------------
  HAlignment ImplAlignatorDots::align(const HAlignandum row, const HAlignandum col, HAlignment result) 
    {
      debug_func_cerr(5);


      startUp(row, col, result);

      performAlignment(row, col, result);

      traceBack(row, col, result);

      cleanUp(row, col, result);

      return result;
    }

  //-----------------------------------------< BackTracke >-------------------------------------------------------------
  void ImplAlignatorDots::traceBack( const HAlignandum row, const HAlignandum col, HAlignment result) 
    {
      debug_func_cerr(5);


      int col_res, row_res; 

      int idot   = mLastDot;
      int jleft  = row->getLength();

      ResiduePair p;

      while ( idot >= 0) {

#ifdef DEBUG
        cout <<
        "-->idot "     << setw(5) << idot      <<
        " col[idot] "  << setw(5) << mPairs[idot].mCol <<
        " row[idot] "  << setw(5) << mPairs[idot].mRow <<
        " mTrace[idot] "<< setw(5) << mTrace[idot] << endl;
#endif

        row_res = mPairs[idot].mRow;
        col_res = mPairs[idot].mCol;

        if (row_res < 1) continue;
        if (col_res < 1) continue;                         
        if (row_res > jleft) break;
        jleft = row_res;                                   // just in case

        result->addPair( p = ResiduePair(row_res, col_res, mPairs[idot].mScore) );

        idot = mTrace[idot];
      }

      result->setScore( mScore );
    } 

  //------------------------------------------------------------------------------------------------------------
  // find the index of a residue pairs given row and column
  Position ImplAlignatorDots::getPairIndex( Position r, Position c ) const {
    int x     = mRowIndices[r];     
    bool found = false;  

    if ( x == NODOT ) 
      return NODOT;  

    while (mPairs[x].mRow == r ) {         
      if (mPairs[x].mCol == c) {             
        found = true;             
        break;         
      }         
      x++;     
    }  

    if (found)         
      return x;     
    else         
      return NODOT; 
  }

  //----------------------------------------------------------------------------------------------------------------------------------------

  Score ImplAlignatorDots::getGapCost( Dot x1, Dot x2 ) const {

    Position c1 = mPairs[x1].mCol;
    Position c2 = mPairs[x2].mCol;  
    Position r1 = mPairs[x1].mRow;
    Position r2 = mPairs[x2].mRow;  

    Score gap_cost = 0;
    Position d;

    if ((d = (r2 - r1)) > 1)
      gap_cost += mRowGop + d * mRowGep;

    if ((d = (c2 - c1)) > 1)
      gap_cost += mColGop + d * mColGep;

    return gap_cost;
  }

  //-----------------------------------------------------------< Alignment subroutine >----------------------------------------------
  void ImplAlignatorDots::performAlignment( const HAlignandum prow, const HAlignandum pcol, HAlignment ali) {

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

    int left;
    int bestdot, globdot;

    Score globbest, best;
    Score sa,sb,sc,sd,s;

    int i;

    int last_search_col, xcol = 0;
    int search_dot, prev_col_dot, prev_row_dot, diagonal_dot, last_dot, current_dot, xdot;

    /* residue numbers of row and column for current dot. These are counted starting at 1 */
    int col_res, row_res;

    int bestpercolstackptr;					/* points to next free element in bestpercolstack */

    /* Allocate Memory for bestpercol and topdot */
    int *bestpercol      = new int[mColLength + 1];
    int *bestpercolstack = new int[mColLength + 1];
    Score *m          = new Score[mNDots];

    glob_dot = search_dot = prev_col_dot = prev_row_dot = diagonal_dot = last_dot = -1;
    globbest = best = 0;

    bestpercolstackptr = STACKEMPTY;

    for (col_res = 1; col_res <= mColLength; col_res++) { bestpercol[col_res] = -1; }
    for (i = 0; i < mNDots; i ++) { mTrace[i]     = -1; m[i] = 0; }

    last_search_col = 1;
    current_dot = 0;

    //#define DEBUG
    //----------------------------------> main alignment loop <----------------------------------------------------
    for ( current_dot = mRowIndices[1]; current_dot < mNDots; current_dot++ ) {	   /* iterate through nextrow starting at first position */

      if (current_dot < 0) continue;
      row_res = mPairs[current_dot].mRow;                           /* row_res = row */
      col_res = mPairs[current_dot].mCol;                           /* col_res = col, wrap around col */

      // some safety checks
#ifdef SAVE
      if (row_res > mRowLength) break;                   /* checking boundaries */
      if (col_res < 1) continue;                /* can be removed, if only dots in boundaries are supplied */
      if (col_res > mColLength) continue;
#endif

#ifdef DEBUG
      printf("--------------------------------------------\n");
      printf("current_dot = %i, row_res = %i, col_res = %i, score = %5.2f\n", current_dot, row_res, col_res, mPairs[current_dot].mScore);
#endif
      /* calculate top row */

      /*------------------------------------------------------------------------------*/
      if ( (last_dot < 0) ||				/* enter first time */
          (row_res  > mPairs[last_dot].mRow) ) {		/* skip, if not in the same row as last time*/

            /* commit changes to bestpercol from bestpercolstack */
            while( bestpercolstackptr > STACKEMPTY ) {
              xdot = bestpercolstack[--bestpercolstackptr];
              xcol = mPairs[xdot].mCol;
              if ( (bestpercol[xcol] == -1) || best > m[bestpercol[xcol]])
                bestpercol[xcol] = xdot;
            }

            /* init all */
            last_search_col = 1; 
            search_dot = -1; prev_col_dot = -1; prev_row_dot = -1; 
      }

      /*------------------------------------------------------------------------------*/
      /* update prev_row_dot = maximum dot along this row */
      xdot = mRowIndices[row_res-1];
      sc = 0; 
      while ( (xdot > -1)  && 
          (mPairs[xdot].mRow == row_res-1 ) &&  /* stop, if dot in previous row any more*/
          (mPairs[xdot].mCol  < col_res-1 )     /* end, if direct contact to new dot*/
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
        last_search_col = col_res;
      } else { 
        search_dot = -1; 
      }

      /*------------------------------------------------------------------------------*/
      /*  select best of d|a|b|c */
      best = 0; bestdot = -1; 

      if( sd > best ) { best = sd; bestdot = diagonal_dot; }
      if( sa > best)  { best = sa; bestdot = search_dot; }
      if( sb > best)  { best = sb; bestdot = prev_col_dot; }
      if( sc > best)  { best = sc; bestdot = prev_row_dot; }

      /* record mTraceback */
      best += mPairs[current_dot].mScore; 

      if (best < 0) { /* local alignment, reset to zero or start new mTrace with single match */
        best    = 0; 
        bestdot = -1;
      }				
      m[current_dot]      = best;
      mTrace[current_dot] = bestdot; 

      if ( best > globbest) { /* save best dot */
        globbest = best; 
        globdot  = current_dot;
      }

#ifdef DEBUG
      printf("current_dot %5i; mTrace[current_dot] %5i; m[current_dot] %5.2f; best %5.2f\n",
          current_dot, mTrace[current_dot], m[current_dot], best);
      printf("search_dot %5i; bdot %5i; prev_row_dot %5i; diagonal_dot %5i\n",
          search_dot,bdot,prev_row_dot,diagonal_dot);
      printf("sa   %5.2f; sb   %5.2f; sc   %5.2f; sd   %5.2f\n",
          sa, sb, sc, sd);
      printf("last_search_col %5i; globdot %5i; globbest %5.2f\n",
          last_search_col, globdot, globbest );
#endif

      /* save score as best score per col, if score is higher than previous score */
      if ( (bestpercol[col_res] == -1) || best > m[bestpercol[col_res]])
        bestpercolstack[bestpercolstackptr++] = current_dot;

      last_dot = current_dot;
    }

    mLastDot= globdot;
    mScore  = globbest;

    //--------------> cleaning up <---------------------------------------------------------------

    delete [] bestpercolstack;
    delete [] bestpercol;        
    delete [] m;        

    //--------------------------------------------------------------------------------------------------
  }

} // namespace alignlib
