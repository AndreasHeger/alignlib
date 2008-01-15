/*
  alignlib - a library for aligning protein sequences

  $Id: ImplDottorSimilarity.cpp,v 1.2 2004/01/07 14:35:35 aheger Exp $

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
#include "alignlib_fwd.h"
#include "alignlib_interfaces.h"
#include "AlignlibDebug.h"
#include "Alignandum.h"
#include "SubstitutionMatrix.h"
#include "ImplDottorSimilarity.h"
#include "ImplAlignmentMatrix.h"

#ifdef WITH_DMALLOC
#include <dmalloc.h>
#endif

using namespace std;

namespace alignlib {

Dottor * makeDottorSimilarity( const SubstitutionMatrix * matrix) {
    return new ImplDottorSimilarity( matrix );
}

//---------------------------------------------------------< constructors and destructors >--------------------------------------
ImplDottorSimilarity::ImplDottorSimilarity ( const SubstitutionMatrix * matrix) : ImplDottor(), mSubstitutionMatrix( matrix ) {
}
		       
ImplDottorSimilarity::~ImplDottorSimilarity () {
}

ImplDottorSimilarity::ImplDottorSimilarity (const ImplDottorSimilarity & src ) : ImplDottor(src), mSubstitutionMatrix( src.mSubstitutionMatrix) {
}

void ImplDottorSimilarity::calculateNewPairs( const HAlignandum row, const HAlignandum col) const {
    Position i, j;

    Position row_length = row->getLength();
    Position col_length = col->getLength();
    
    ResiduePAIR p;

    for (i = 1; i <= row_length; i++) 
      for (j = 1; j <= col_length; i++) 
	  if (mSubstitutionMatrix->getScore(row->asResidue(i), row->asResidue(j)) > 0 )
	      mMatrix->addPair( p = ResiduePAIR( i, j, 1 ));
    
};


} // namespace alignlib
