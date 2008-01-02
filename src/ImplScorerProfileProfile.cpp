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
#include "ImplProfile.h"
#include "ImplScorerProfileProfile.h"
#include "ImplTranslator.h" // for direct access to mask_code

using namespace std;

namespace alignlib
{

  // factory function for creating iterator over full matrix
  Scorer * makeScorerProfileProfile( const Alignandum * row, const Alignandum * col )
  {
    return new ImplScorerProfileProfile( row, col );
  }
  
  //--------------------------------------------------------------------------------------
  ImplScorerProfileProfile::ImplScorerProfileProfile( const Alignandum * row,
						      const Alignandum * col ) :
    ImplScorer( row, col )
  {
    const ImplProfile * s1 = dynamic_cast<const ImplProfile*>(row);
    const ImplProfile * s2 = dynamic_cast<const ImplProfile*>(col);	

    if (!s1)
      throw AlignException( "ImplScoreProfileProfile.cpp: row not a profile.");
    
    if (!s2)
      throw AlignException( "ImplScoreProfileProfile.cpp: col not a profile.");
    
    mRowProfile     = s1->getData().mProfilePointer;
    mRowFrequencies = s1->getData().mFrequenciesPointer;    
    
    mColProfile     = s2->getData().mProfilePointer;
    mColFrequencies = s2->getData().mFrequenciesPointer;        

    if ( s1->getTranslator()->getAlphabetSize() != 
    	s2->getTranslator()->getAlphabetSize() )
    	throw AlignException( "ImplScorerProfileProfile.cpp: alphabet size different in row and col");
        
    mAlphabetSize = s1->getTranslator()->getAlphabetSize();
  }    
  
  //--------------------------------------------------------------------------------------
  ImplScorerProfileProfile::~ImplScorerProfileProfile ()
  {
  }
  
  //--------------------------------------------------------------------------------------
  ImplScorerProfileProfile::ImplScorerProfileProfile(const ImplScorerProfileProfile & src) :
    ImplScorer(src),
    mRowProfile( src.mRowProfile ), mRowFrequencies( src.mRowFrequencies ),
    mColProfile( src.mColProfile ), mColFrequencies( src.mColFrequencies )
  {
  }

  //--------------------------------------------------------------------------------------  
  /** return a copy of the same iterator
   */
  Scorer * ImplScorerProfileProfile::getClone() const
  {
    return new ImplScorerProfileProfile( *this );
  }

  /** return a new instance of this object initialized with row and col
   */
  Scorer * ImplScorerProfileProfile::getNew( const Alignandum * row, const Alignandum * col) const
  {
    return new ImplScorerProfileProfile( row, col );
  }
  
  /** return score of matching row to col
   */
  Score ImplScorerProfileProfile::getScore( Position row, Position col ) const
  {
    Score score = 0;
    for (int i = 0; i < mAlphabetSize; i++ ) 
      score +=	mRowProfile[row][i] * mColFrequencies[col][i] + 
	mColProfile[col][i] * mRowFrequencies[row][i];
    
    return score;
  }
  
}

  
