/*
  alignlib - a library for aligning protein sequences

  $Id: Iterator2D.cpp,v 1.2 2004/01/07 14:35:32 aheger Exp $

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

#include "alignlib_fwd.h"
#include "alignlib_interfaces.h"
#include "alignlib_fwd.h"
#include "AlignException.h"
#include "Alignandum.h"
#include "Translator.h"
#include "ImplSequence.h"
#include "ImplProfile.h"
#include "ImplScorerSequenceProfile.h"

using namespace std;

namespace alignlib
{

  // factory function for creating iterator over full matrix
  HScorer makeScorerSequenceProfile( 
		  const HAlignandum & row, 
		  const HAlignandum & col )
  {
    return HScorer( new ImplScorerSequenceProfile( row, col ) );
  }
  
  //--------------------------------------------------------------------------------------
  ImplScorerSequenceProfile::ImplScorerSequenceProfile( 
		  const HAlignandum & row,
		  const HAlignandum & col ) :
    ImplScorer( row, col )
  {
    const HImplSequence s1 = boost::dynamic_pointer_cast<ImplSequence, Alignandum>(row);	
    const HImplProfile s2 = boost::dynamic_pointer_cast<ImplProfile, Alignandum>(col);

    mRowSequence    = s1->getSequence();
    mColProfile     = s2->getScoreMatrix();

    mProfileWidth = s2->getTranslator()->getAlphabetSize();
    
    if ( s1->getTranslator()->getAlphabetSize() != mProfileWidth ) 
    	throw AlignException( "ImplScorerSequenceProfile.cpp: alphabet size different in row and col");
    
  }    
  
  //--------------------------------------------------------------------------------------
  ImplScorerSequenceProfile::~ImplScorerSequenceProfile ()
  {
  }
  
  //--------------------------------------------------------------------------------------
  ImplScorerSequenceProfile::ImplScorerSequenceProfile(const ImplScorerSequenceProfile & src) :
    ImplScorer(src),
    mRowSequence( src.mRowSequence ),
    mColProfile( src.mColProfile ),
    mProfileWidth( src.mProfileWidth )
  {
  }

  //--------------------------------------------------------------------------------------  
  /** return a copy of the same iterator
   */
  HScorer ImplScorerSequenceProfile::getClone() const
  {
    return HScorer( new ImplScorerSequenceProfile( *this ) );
  }

  /** return a new instance of this object initialized with row and col
   */
  HScorer ImplScorerSequenceProfile::getNew( 
		  const HAlignandum & row, 
		  const HAlignandum & col) const
  {
    return HScorer( new ImplScorerSequenceProfile( row, col ) );
  }

  /** return score of matching row to col
   */
  Score ImplScorerSequenceProfile::getScore( Position row, Position col ) const
  {
	return mColProfile->getValue( col, mRowSequence[row] );
  }
  
}

  
