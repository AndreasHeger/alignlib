/*
  alignlib - a library for aligning protein sequences

  $Id: ImplAlignatorFullDPWrap.cpp,v 1.2 2004/01/07 14:35:34 aheger Exp $

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
#include "ImplAlignatorFullDPWrap.h"
#include "HelpersAlignator.h"
#include "HelpersSubstitutionMatrix.h"

#ifdef WITH_DMALLOC
#include <dmalloc.h>
#endif

using namespace std;

namespace alignlib {

//--------------------------------> factory functions <--------------------------------------------------------------------
Alignator * makeAlignatorFullDPWrap( Score gop, Score gep, const SubstitutionMatrix * matrix ) {
  if (matrix) 
    return new ImplAlignatorFullDPWrap( matrix, gop, gep, gop, gep);
  else
    return new ImplAlignatorFullDPWrap( getDefaultSubstitutionMatrix(), gop, gep, gop, gep);
}

//----------------------------------------------------------------------------------------------------------------------------------------
ImplAlignatorFullDPWrap::ImplAlignatorFullDPWrap( const SubstitutionMatrix * subst_matrix,
					  Score row_gop, Score row_gep, 
					  Score col_gop, Score col_gep ) :
    ImplAlignatorFullDP( subst_matrix, row_gop, row_gep, col_gop, col_gep) {
}
  
//----------------------------------------------------------------------------------------------------------------------------------------
ImplAlignatorFullDPWrap::ImplAlignatorFullDPWrap( const ImplAlignatorFullDPWrap & src ) : ImplAlignatorFullDP( src ) 
{
  debug_func_cerr(5);

}

//----------------------------------------------------------------------------------------------------------------------------------------
ImplAlignatorFullDPWrap::~ImplAlignatorFullDPWrap() 
{
  debug_func_cerr(5);


}

//----------------------------------------------------------------------------------------------------------------------------------------
void ImplAlignatorFullDPWrap::startUp(const Alignandum * row, const Alignandum *col, Alignata * ali) {
    ImplAlignatorFullDP::startUp(row, col, ali);  
    
    for (int r = 0; r <= mRowLength; r++)
	mTrace[getTraceIndex(r,0)] = TB_WRAP;
}

//---------------------------------< the actual alignment algorithm >-------------------------------------------
void ImplAlignatorFullDPWrap::performAlignment( const Alignandum * prow, const Alignandum * pcol, Alignata * ali) {

  register int   row, col;

  Position row_length = getRowLength();
  Position col_length = getColLength();

  Score row_gop = getRowGop();
  Score row_gep = getRowGep();
  Score col_gop = getColGop();
  Score col_gep = getColGep();

  mRowLast = 0;
  mColLast = 0;
  mScore = 0;

  Score c, e, d, s;                  // helper variables
  
  Score row_m = row_gop + row_gep;
  Score col_m = col_gop + col_gep;

  int matchtype;

  //----------------------------> Initialise affine penalty arrays <-------------------------------
    mCC[0] = 0;
    for (col = 1; col <= col_length; col++) {
      mCC[col]   = 0;
      mDD[col]   = row_gop;                               // score for horizontal gap opening
    }
 
    mCC[col_length]   = col_gop;
    
    //----------------------------> Calculate dynamic programming matrix <----------------------------
    //----------------------------> iterate over rowumns <--------------------------------------------
    for (row = 1; row <= row_length; row++) {

	// this part is different from the ordinary alignment
      matchtype = TB_MATCH;

      if (mCC[col_length] > 0) {					// the wrapping around part
	mCC[0] = c = mCC[col_length];
      } else {
	mCC[0] = c = 0;
      }

      s = mCC[0];
      
      e = col_gop;                                        // penalty for opening a vertical gap
	
      //-------------------------> iterate over cols <------------------------------------------------
      for (col = 1; col <= col_length; col++) {
 
	    //---------------------------> calculate scores <--------------------------------------------
	    // c contains score of cell above
	    // s contains score for cell [row, col-1]
	    // e is better of: score for opening a vertical gap or score for extending a vertical gap: use col-gap-penalties
	    if ((c =   c     + col_m) > (e =     e   + col_gep))  e = c;
	    // d is better of: score for opening a horizontal gap or score for extending a horizontal gap
	    if ((c = mCC[col] + row_m) > (d = mDD[col] + row_gep))  d = c;
	    
	    // c is score for a match
	    c = s + (this->*mMatchFunction)( row, col );
	    // put into c the best of all possible cases
	    if (e > c) c = e;
	    if (d > c) c = d;
	    
	    //--------------------------> recurrence relation <-------------------------------------------------
	    if (c <= 0) {
	      c = 0;                                                  // the local alignment part
	    } else {
	      if ( c == d )                   // horizontal gap
		mTrace[getTraceIndex(row,col)] = TB_INSERTION;
	      else if ( c == e )              // vertical gap
		mTrace[getTraceIndex(row,col)] = TB_DELETION;
	      else {                          // match
		mTrace[getTraceIndex(row,col)] = matchtype;
		matchtype = TB_MATCH;
	      }
	    }

	    s = mCC[col];
	    mCC[col] = c;                                              // save new score for next i
	    mDD[col] = d;
	    
	    if (mScore < c) {                                            // save maximum
	      mScore   = c;
	      mRowLast = row;
	      mColLast = col;
	    }
      }
    }

}   



} // namespace alignlib
