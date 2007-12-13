/*
  alignlib - a library for aligning protein sequences

  $Id: HelpersAlignmentPrivate.cpp,v 1.3 2004/06/02 12:11:37 aheger Exp $

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
#include <string>
#include <stdio.h>
#include "alignlib.h"
#include "AlignlibDebug.h"

#include "Alignandum.h"
#include "AlignmentIterator.h"
#include "Alignment.h"
#include "AlignException.h"

#include "HelpersAlignment.h"
#include "Translator.h"
#include "HelpersTranslator.h"
#include "ImplTranslator.h"

#include "Alignatum.h"
#include "HelpersAlignatum.h"

#include "SubstitutionMatrix.h"
#include "HelpersSubstitutionMatrix.h"

#include "ImplSequence.h"
#include "ImplProfile.h"

#include "AlignException.h"

#ifdef WITH_DMALLOC
#include <dmalloc.h>
#endif

using namespace std;

namespace alignlib {

  // a class that calculates the match-score between a row and column.

    class Match {

	typedef Score (Match::*MATCH_FUNCTION_POINTER)( Position, Position); 
    public:
	
	Match( const Alignandum * row, const Alignandum * col, const SubstitutionMatrix * matrix = NULL) {
	    
	    if (matrix)
		mSubstitutionMatrix = matrix;
	    else 
		mSubstitutionMatrix = getDefaultSubstitutionMatrix();
	    
	    //-------------------------------------------------------------------------------------------
	    // try to typecast row and col
	    // const SequenceTCO & t1 = dynamic_cast<const SequenceTCO&>(row);
	    const ImplSequence * s1 = dynamic_cast<const ImplSequence*>(row);
	    const ImplProfile * p1 = dynamic_cast<const ImplProfile*>(row);
	    
	    // const SequenceTCO & t2 = dynamic_cast<const SequenceTCO&>(col);
	    const ImplSequence * s2 = dynamic_cast<const ImplSequence*>(col);
	    const ImplProfile * p2 = dynamic_cast<const ImplProfile*>(col);
	    
	    // setup static pointers to the data locations
	    if (s1) {
		mRowSequence = s1->getData().mSequencePointer;
	    } else if (p1) {
	        mRowProfile     = p1->getData().mProfilePointer;
		mRowFrequencies = p1->getData().mFrequenciesPointer;
	    }
	    
	    if (s2) {
		mColSequence = s2->getData().mSequencePointer;
	    } else if (p2) {
		mColProfile     = p2->getData().mProfilePointer;
		mColFrequencies = p2->getData().mFrequenciesPointer;
	    }
	    
	    //!! to do: add TCO
	    if (s1 && s2) {
		mMatchFunction = &Match::matchSequenceSequence;
	    }
	    
	    if (p1 && p2) {
		mMatchFunction = &Match::matchProfileProfile;
	    }

	    if (s1 && p2) {
	      mMatchFunction = &Match::matchSequenceProfile;
	    }

	    if (s2 && p1) {
	      mMatchFunction = &Match::matchProfileSequence;
	    }
	    
	    // throw exception if you can not handle row/col combination
	    //!! to be implemented
	    
	    if (!row->isPrepared() )
		throw AlignException( "Row has not been prepared in HelpersAlignmentPrivate::Match!");

	    if (!col->isPrepared() ) 
		throw AlignException( "Col has not been prepared in HelpersAlignmentPrivate::Match!");
	}

	/** call matchfunction for row and col */
	Score operator()(Position row, Position col) {
	    return (this->*mMatchFunction)( row, col );
	}

    private:
	//------------------------------------< match functions >--------------------------------------------
	Score matchSequenceSequence( Position row, Position col ) {
	    return (mSubstitutionMatrix->getScore(mRowSequence[row],mColSequence[col]));
	}
	
	Score matchProfileSequence( Position row, Position col ) {
	  if (mColSequence[col] == CODE_MASK)
	    return MASK_VALUE;
	  else
	    return mRowProfile[row][mColSequence[col]];
	}
	
	Score matchSequenceProfile( Position row, Position col ) {
	  if (mRowSequence[row] == CODE_MASK)
	    return MASK_VALUE;
	  else
	    return mColProfile[col][mRowSequence[row]];
	}
	
	Score matchProfileProfile( Position row, Position col ) {
	    Score score = 0;
	    for (Position i = 0; i < PROFILEWIDTH; i++ ) 
		score += mRowProfile[row][i] * mColFrequencies[col][i] + 
		         mColProfile[col][i] * mRowFrequencies[row][i];
	    return score;
	}

	/** pointer to function that calculates match-score */
	MATCH_FUNCTION_POINTER mMatchFunction;
	
	/** pointer to member data of row/col : AlignandumSequence */
	const Residue * mRowSequence; 
	const Residue * mColSequence;
	
	/** pointer to member data of row/col AlignandumProfile */
	const ProfileColumn * mRowProfile;
	const ProfileColumn * mColProfile;
	const FrequencyColumn * mRowFrequencies;
	const FrequencyColumn * mColFrequencies;
	
	/** pointer to substitution matrix to use */
	const SubstitutionMatrix * mSubstitutionMatrix;

    };


//----------------------------------------------------------------------------------------------------
Alignment * rescoreAlignmentPrivate( Alignment * dest,
				    const Alignandum * row,
				    const Alignandum * col, 
				    const SubstitutionMatrix * matrix) 
{
  debug_func_cerr(5);

  
  Match match( row, col, matrix);

  AlignmentIterator it(dest->begin());
  AlignmentIterator it_end(dest->end());

  for (; it != it_end; ++it) 
      it->mScore = match( it->mRow, it->mCol );

  return dest;
}


} // namespace alignlib
