/*
  alignlib - a library for aligning protein sequences

  $Id: ImplFragmentorRepetitive.cpp,v 1.2 2004/01/07 14:35:35 aheger Exp $

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
#include "alignlib_fwd.h"
#include "AlignlibDebug.h"

#include "Fragmentor.h"

#include "Alignment.h"
#include "HelpersAlignment.h"

#include "AlignmentIterator.h"

#include "Alignandum.h"
#include "AlignException.h"

#include "HelpersSubstitutionMatrix.h"

#include "Alignator.h"
#include "HelpersAlignator.h"

#include "ImplFragmentorRepetitive.h"

using namespace std;

namespace alignlib 
{

/*---------------------factory functions ---------------------------------- */

  /** make an alignator object, which does a dot-alignment. The default version can be given an AlignmentMatrix-
      object */
Fragmentor * makeFragmentorRepetitive( Alignator * alignator, Score min_score ) {
  return new ImplFragmentorRepetitive( alignator, min_score );
}


//----------------------------------------------------------------------------------------------------
  
ImplFragmentorRepetitive::ImplFragmentorRepetitive( Alignator * alignator, 
						    Score min_score ):

  ImplFragmentor(),
  mAlignator( alignator ),
  mMinScore( min_score) {
}


ImplFragmentorRepetitive::~ImplFragmentorRepetitive() 
{
  debug_func_cerr(5);

}

ImplFragmentorRepetitive::ImplFragmentorRepetitive( const ImplFragmentorRepetitive & src ) : 
    ImplFragmentor(src), 
    mAlignator(src.mAlignator), 
    mMinScore( src.mMinScore ) {
}

//------------------------------------------------------------------------------------------------
void ImplFragmentorRepetitive::performFragmentation( const Alignandum * row, 
						     const Alignandum * col, 
						     const Alignment * sample) {
    
  /* since src1 and src2 are const, I have to create two work-copies, 
     so that the boundaries can be changed. */

  Alignandum * copy_row = row->getClone();
  Alignandum * copy_col = col->getClone();
  
  while ( 1 ) {
    
    Alignment * result = sample->getNew();
    mAlignator->align( copy_row, copy_col, result);
    if (result->getScore() >= mMinScore) {

      mFragments->push_back( result );    
      copy_row->mask( result->getRowFrom(), result->getRowTo() );
      copy_col->mask( result->getColFrom(), result->getColTo() );
      
    } else {
      delete result;    
      break;      
    }
  }

  delete copy_row;
  delete copy_col;

}
  

} // namespace alignlib




