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
#include "AlignException.h"
#include "Alignandum.h"
#include "Translator.h"
#include "ImplSequence.h"
#include "ImplProfile.h"
#include "ImplScorerProfileSequence.h"

using namespace std;

namespace alignlib
{

  // factory function for creating iterator over full matrix
  Scorer * makeScorerProfileSequence( const Alignandum * row, const Alignandum * col )
  {
    return new ImplScorerProfileSequence( row, col );
  }
  
  //--------------------------------------------------------------------------------------
  ImplScorerProfileSequence::ImplScorerProfileSequence( const Alignandum * row,
							const Alignandum * col ) :
    ImplScorer( row, col )
  {
    const ImplProfile * s1 = dynamic_cast<const ImplProfile*>(row);
    const ImplSequence * s2 = dynamic_cast<const ImplSequence*>(col);	

    if (!s1)
      throw AlignException( "ImplScoreProfileSequence.cpp: row not a profile.");
    
    if (!s2)
      throw AlignException( "ImplScoreProfileSequence.cpp: col not a sequence.");
        
    mRowProfile     = s1->getData().mProfilePointer;
    mColSequence    = s2->getData().mSequencePointer;

    if ( s1->getTranslator()->getAlphabetSize() != 
    	s2->getTranslator()->getAlphabetSize() )
    	throw AlignException( "ImplScorerProfileSequence.cpp: alphabet size different in row and col");
    
  }    
  
  //--------------------------------------------------------------------------------------
  ImplScorerProfileSequence::~ImplScorerProfileSequence ()
  {
  }
  
  //--------------------------------------------------------------------------------------
  ImplScorerProfileSequence::ImplScorerProfileSequence(const ImplScorerProfileSequence & src) :
    ImplScorer(src),
    mRowProfile( src.mRowProfile ),
    mColSequence( src.mColSequence )
  {
  }

  //--------------------------------------------------------------------------------------  
  /** return a copy of the same iterator
   */
  Scorer * ImplScorerProfileSequence::getClone() const
  {
    return new ImplScorerProfileSequence( *this );
  }

  /** return a new instance of this object initialized with row and col
   */
  Scorer * ImplScorerProfileSequence::getNew( const Alignandum * row, const Alignandum * col) const
  {
    return new ImplScorerProfileSequence( row, col );
  }

  /** return score of matching row to col
   */
  Score ImplScorerProfileSequence::getScore( Position row, Position col ) const
  {
	  return mRowProfile[row][mColSequence[col]];
  }
  
}

  
