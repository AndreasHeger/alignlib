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
#include "ImplScorerSequenceSequence.h"
#include "SubstitutionMatrix.h"
#include "ImplSequence.h"

#ifdef WITH_DMALLOC
#include <dmalloc.h>
#endif

using namespace std;

namespace alignlib
{

  // factory function for creating iterator over full matrix
  Scorer * makeScorerSequenceSequence( const Alignandum * row, const Alignandum * col, const SubstitutionMatrix * matrix )
  {
    return new ImplScorerSequenceSequence( row, col, matrix );
  }
  
  //--------------------------------------------------------------------------------------
  ImplScorerSequenceSequence::ImplScorerSequenceSequence( const Alignandum * row,
							  const Alignandum * col,
							  const SubstitutionMatrix * matrix) :
    ImplScorer( row, col )
  {
    const ImplSequence * s1 = dynamic_cast<const ImplSequence*>(row);	
    const ImplSequence * s2 = dynamic_cast<const ImplSequence*>(col);

    if (!s1)
      throw AlignException( "ImplScoreSequenceSequence.cpp: row not a sequence.");
    
    if (!s2)
      throw AlignException( "ImplScoreSequenceSequence.cpp: col not a sequence.");
    
    mRowSequence = s1->getData().mSequencePointer;
    mColSequence = s2->getData().mSequencePointer;

    mSubstitutionMatrix = matrix;
  }    
  
  //--------------------------------------------------------------------------------------
  ImplScorerSequenceSequence::~ImplScorerSequenceSequence ()
  {
  }
  
  //--------------------------------------------------------------------------------------
  ImplScorerSequenceSequence::ImplScorerSequenceSequence(const ImplScorerSequenceSequence & src) :
    ImplScorer(src ),
    mRowSequence( src.mRowSequence ), mColSequence( src.mColSequence ),
    mSubstitutionMatrix( src.mSubstitutionMatrix )
  {
  }

  //--------------------------------------------------------------------------------------  
  /** return a copy of the same iterator
   */
  Scorer * ImplScorerSequenceSequence::getClone() const
  {
    return new ImplScorerSequenceSequence( *this );
  }

  /** return a new instance of this object initialized with row and col
   */
  Scorer * ImplScorerSequenceSequence::getNew( const Alignandum * row, const Alignandum * col) const
  {
    return new ImplScorerSequenceSequence( row, col, mSubstitutionMatrix );
  }
  
  /** return score of matching row to col
   */
  Score ImplScorerSequenceSequence::getScore( Position row, Position col ) const
  {
    return mSubstitutionMatrix->getScore(mRowSequence[row],mColSequence[col]);
  }
  
}

  
