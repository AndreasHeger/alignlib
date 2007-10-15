/*
  alignlib - a library for aligning protein sequences

  $Id: ImplAlignatorFullDPGlobal.cpp,v 1.4 2005/02/24 11:07:25 aheger Exp $

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
#include "AlignException.h"
#include "ImplAlignatorFullDPGlobal.h"
#include "HelpersAlignator.h"
#include "HelpersSubstitutionMatrix.h"

#ifdef WITH_DMALLOC
#include <dmalloc.h>
#endif

using namespace std;

namespace alignlib {

//--------------------------------> factory functions <--------------------------------------------------------------------
Alignator * makeAlignatorFullDPGlobal( Score gop, Score gep, bool penalize_left, bool penalize_right, const SubstitutionMatrix * matrix ) {
  if (matrix) 
    return new ImplAlignatorFullDPGlobal( matrix, gop, gep, gop, gep, penalize_left, penalize_right, penalize_left, penalize_right);
  else
    return new ImplAlignatorFullDPGlobal( getDefaultSubstitutionMatrix(), gop, gep, gop, gep, penalize_left, penalize_right, penalize_left, penalize_right);
}

//----------------------------------------------------------------------------------------------------------------------------------------
ImplAlignatorFullDPGlobal::ImplAlignatorFullDPGlobal( const SubstitutionMatrix * subst_matrix,
						      Score row_gop, Score row_gep, 
						      Score col_gop, Score col_gep,
						      bool penalize_row_left,
						      bool penalize_row_right,
						      bool penalize_col_left,
						      bool penalize_col_right ) :
  ImplAlignatorFullDP( subst_matrix, row_gop, row_gep, col_gop, col_gep),
  mPenalizeRowLeft( penalize_row_left ),
  mPenalizeRowRight( penalize_row_right ),
  mPenalizeColLeft( penalize_col_left ),
  mPenalizeColRight( penalize_col_right ) {
}
  
//----------------------------------------------------------------------------------------------------------------------------------------
ImplAlignatorFullDPGlobal::ImplAlignatorFullDPGlobal( const ImplAlignatorFullDPGlobal & src ) : 
  ImplAlignatorFullDP( src ),
  mPenalizeRowLeft( src.mPenalizeRowLeft ),
  mPenalizeRowRight( src.mPenalizeRowRight ),
  mPenalizeColLeft( src.mPenalizeColLeft ),
  mPenalizeColRight( src.mPenalizeColRight ) 
{
  debug_func_cerr(5);

}

//----------------------------------------------------------------------------------------------------------------------------------------
ImplAlignatorFullDPGlobal::~ImplAlignatorFullDPGlobal() 
{
  debug_func_cerr(5);


}

//---------------------------------< the actual alignment algorithm >-------------------------------------------
void ImplAlignatorFullDPGlobal::performAlignment( const Alignandum * prow, const Alignandum * pcol, Alignata * ali) {

  register int   row, col;

  Position row_length = getRowLength();
  Position col_length = getColLength();

  Score row_gop = getRowGop();
  Score row_gep = getRowGep();
  Score col_gop = getColGop();
  Score col_gep = getColGep();

  mRowLast = 0;
  mColLast = 0;
  mScore = -9999;

  Score c, e, d, s;                  // helper variables
  
  Score row_m = row_gop + row_gep;
  Score col_m = col_gop + col_gep;

  /*
   ------> col
   | CC-> 
   | DD->
   |
   |
   row

   For each cell:
   s     mCC/mDD   
      \  |
       \ |
   c/e-- x
   
   c/mCC: last op was match
   e/mDD: last op was gap

  */
  
  //----------------------------> Initialise affine penalty arrays <-------------------------------
  /* set initial values for upper border */
  c = 0;

  if (mPenalizeRowLeft) {
    for (col = 1; col <= col_length; col++) {
      mCC[col]   = row_gop + row_gep * col;
      mDD[col]   = mCC[col];                              // add score for gap opening
    }
  } else {
    for (col = 1; col <= col_length; col++) {
      mCC[col]   = 0;
      mDD[col]   = row_gop;                               // add score for gap opening
    }
  }    

  
  //----------------------------> Calculate dynamic programming matrix <----------------------------
  //----------------------------> iterate over rowumns <--------------------------------------------
  for (row = 1; row <= row_length; row++) {

    /* set initial values for left border */
    if (mPenalizeColLeft) {
      e = c = col_gop + col_gep * row;				
      s = (row - 1) * col_gep;
      if (row > 1) 
	s += row_gop;

    } else {
      s = 0;
      c = 0;
      e = col_gop;
    }
    
    //-------------------------> iterate over cols <------------------------------------------------
    for (col = 1; col <= col_length; col++) {
      
      // cout << "row=" << row << " col=" << col << " c=s= " << s << " e=" << e << endl;
      
      //---------------------------> calculate scores <--------------------------------------------
      // e is better of: score for opening a vertical gap or score for extending a vertical gap: use col-gap-penalties
      if ((c =   c     + col_m) > (e =     e   + col_gep))  e = c;
      // d is better of: score for opening a horizontal gap or score for extending a horizontal gap
      if ((c = mCC[col] + row_m) > (d = mDD[col] + row_gep))  d = c;
      
      // c is score for a match
      c = s + (this->*mMatchFunction)( row, col );
      // put into c the best of all possible cases
      if (e > c) c = e;
      if (d > c) c = d;
      
      // cout << "c=" << c << " d=" << d << " e=" << e << endl;
      
      if ( c == d )                   // horizontal gap
	mTrace[getTraceIndex(row,col)] = TB_INSERTION;
      else if ( c == e )              // vertical gap
	mTrace[getTraceIndex(row,col)] = TB_DELETION;
      else                            // match
	mTrace[getTraceIndex(row,col)] = TB_MATCH;
      
      s = mCC[col];
      mCC[col] = c;                                              // save new score for next i
      mDD[col] = d;

      /* retrieve maximum score. If penalty has to be paid for right end-gaps, 
	 then restrict search correspondingly;
      */
      if (mPenalizeRowRight && row < row_length) 
	continue;
      if (mPenalizeColRight && col < col_length) 
	continue;
      
      if (mScore < c) {                                            // save maximum
	mScore   = c;
	mRowLast = row;
	mColLast = col;
      }
    }
  }
}   


} // namespace alignlib
