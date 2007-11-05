/*
  alignlib - a library for aligning protein sequences

  $Id: HelpersAlignator.cpp,v 1.4 2005/02/24 11:07:25 aheger Exp $

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
#include <math.h>
#include "alignlib.h"
#include "AlignlibDebug.h"
#include "AlignlibDebug.h"
#include "Alignandum.h"

#include "Alignata.h"
#include "HelpersAlignata.h"

#include "HelpersSubstitutionMatrix.h"
#include "HelpersAlignator.h"
#include "ImplAlignatorSimilarity.h"
#include "ImplAlignatorIdentity.h"
#include "ImplAlignatorDPFull.h"


#ifdef WITH_DMALLOC
#include <dmalloc.h>
#endif

using namespace std;

namespace alignlib {

/*---------------------factory functions ---------------------------------- */
Alignator * makeAlignatorIdentity() {
  return new ImplAlignatorIdentity();
}

  /*
Alignator * makeFullDP( Score gop, Score gep, const SubstitutionMatrix * matrix ) {
  return makeAlignatorFullDP( gop, gep, matrix );
}
  
Alignator * makeAlignatorFullDP( Score gop, Score gep, const SubstitutionMatrix * matrix ) {
  if (matrix) 
    return new ImplAlignatorFullDP( matrix, gop, gep, gop, gep);
  else
    return new ImplAlignatorFullDP( getDefaultSubstitutionMatrix(), gop, gep, gop, gep);
}
  */
  
//---------------------------------------------------------------------------
// convenience functions:

Alignata * performIterativeAlignmentStep( Alignata * dest, 
				      Alignandum * copy_1, 
				      Alignandum * copy_2, 
				      Alignator * alignator, 
				      Score min_score) {
  debug_func_cerr(5);
  
  // do alignment with current boundaries of the objects, but remember them
  Position from_1 = copy_1->getFrom();
  Position from_2 = copy_2->getFrom();

  Position to_1 = copy_1->getTo();
  Position to_2 = copy_2->getTo();

  debug_cerr(5, "aligning in regions (" << from_1 << "-" << to_1 <<") -> (" << from_2 << "-" << to_2 << ")" );

  if (from_1 > to_1 || from_2 > to_2)
    return dest;

  Alignata * result = dest->getNew();

  alignator->align( copy_1, copy_2, result );

  if (result->getScore() > min_score) {

      addAlignata2Alignata( dest, result );

      debug_cerr( 5, "new alignment\n" << *result )
      debug_cerr( 5, "new alignment coordinates: row=" << result->getRowFrom() << " " << result->getRowTo()  
    		  	<< " col=" << result->getColFrom() << " " << result->getColTo() );      
      debug_cerr( 5, "current alignment\n" << *result )      

      Position from_1_result = result->getRowFrom();
      Position from_2_result = result->getColFrom();
      Position to_1_result   = result->getRowTo();
      Position to_2_result   = result->getColTo();
      
      // align in region before current alignment
      copy_1->useSegment( from_1, from_1_result);
      copy_2->useSegment( from_2, from_2_result);
      performIterativeAlignmentStep( dest, copy_1, copy_2, alignator, min_score);
      
      // align in region after current alignment
      copy_1->useSegment( to_1_result, to_1);
      copy_2->useSegment( to_2_result, to_2);
      performIterativeAlignmentStep( dest, copy_1, copy_2, alignator, min_score);

  }
  delete result;

  return dest;
} 

// iteratively align src1 and src2, until score drops below min_score
Alignata * performIterativeAlignment( Alignata * dest, 
				      const Alignandum * src_1, 
				      const Alignandum * src_2, 
				      Alignator * alignator, 
				      Score min_score) 
{
  debug_func_cerr(5);


    /* since src1 and src2 are const, I have to create two work-copies, 
       so that the boundaries can be changed. */

  Alignandum * copy_1 = src_1->getClone();
  Alignandum * copy_2 = src_2->getClone();
  
  // start aligning by calling recursively performIterativeAlignmentStep
  performIterativeAlignmentStep( dest, copy_1, copy_2, alignator, min_score );

  delete copy_1;
  delete copy_2;

  return dest;
} 

// iteratively align src1 and src2, until score drops below min_score
Alignata * performIterativeAlignment( Alignata * dest, 
				      Alignandum * src_1, 
				      Alignandum * src_2, 
				      Alignator * alignator, 
				      Score min_score) 
{
  debug_func_cerr(5);


  // do alignment with current boundaries of the objects, but remember them
  Position from_1 = src_1->getFrom();
  Position from_2 = src_2->getFrom();

  Position to_1 = src_1->getTo();
  Position to_2 = src_2->getTo();
  
  // start aligning by calling recursively performIterativeAlignmentStep
  performIterativeAlignmentStep( dest, src_1, src_2, alignator, min_score );

  // restore old segments
  src_1->useSegment( from_1, to_1);
  src_2->useSegment( from_2, to_2);

  return dest;
} 


} // namespace alignlib
