/*
  alignlib - a library for aligning protein sequences

  $Id: ImplDottor.cpp,v 1.2 2004/01/07 14:35:35 aheger Exp $

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
#include "ImplDottor.h"
#include "ImplAlignmentMatrix.h"

#ifdef WITH_DMALLOC
#include <dmalloc.h>
#endif

using namespace std;

namespace alignlib {

    static const Dottor * DEFAULT_DOTTOR = makeDottorIdentity();

    /** gets the default Dottor object */ 
    const Dottor * getDefaultDottor() {
      return DEFAULT_DOTTOR;
    }

  /** sets the default Dottor object */
  const Dottor * setDefaultDottor( const Dottor * dottor ) {
    const Dottor * t = DEFAULT_DOTTOR;
    DEFAULT_DOTTOR = dottor;
    return t;
  }


//---------------------------------------------------------< constructors and destructors >--------------------------------------
ImplDottor::ImplDottor () {
}
		       
ImplDottor::~ImplDottor () {
}

ImplDottor::ImplDottor (const ImplDottor & src ) {
}

/** release memory hold by ImplDottor */
void ImplDottor::releasePairs() const {
  delete mMatrix;
};
    
/** calculate the dots */
void ImplDottor::calculatePairs( const HAlignandum row, const HAlignandum col) const {

  mMatrix = new ImplAlignmentMatrix();

  calculateNewPairs( row, col );

  mMatrix->calculateLength();

};
    
/** get array of residue pairs */
ResiduePair * ImplDottor::getPairs() const {
  return mMatrix->mPairs;
}

/** get location of array for row-indices */
Position * ImplDottor::getRowIndices() const {
  return mMatrix.mRowIndices;
}

} // namespace alignlib
