/*
  alignlib - a library for aligning protein sequences

  $Id: ImplFragmentorIterative.cpp,v 1.2 2004/01/07 14:35:35 aheger Exp $

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
#include "AlignlibDebug.h"

#include "Fragmentor.h"

#include "Alignata.h"
#include "HelpersAlignata.h"

#include "AlignataIterator.h"

#include "Alignandum.h"
#include "AlignException.h"

#include "SubstitutionMatrix.h"
#include "HelpersSubstitutionMatrix.h"

#include "Alignator.h"
#include "HelpersAlignator.h"

#include "ImplFragmentorIterative.h"

#ifdef WITH_DMALLOC
#include <dmalloc.h>
#endif

using namespace std;

namespace alignlib {

/*---------------------factory functions ---------------------------------- */

  /** make an alignator object, which does a dot-alignment. The default version can be given an AlignataMatrix-
      object */
Fragmentor * makeFragmentorIterative( Alignata * dots, Score min_score, Score gop, Score gep) {
  
  return new ImplFragmentorIterative( dots, min_score, gop, gep);
}


//----------------------------------------------------------------------------------------------------
  
ImplFragmentorIterative::ImplFragmentorIterative( Alignata * dots, 
						  Score min_score, 
						  Score gop, 
						  Score gep) :
  ImplFragmentor(),
  mDots(dots),
  mMinScore( min_score),
  mGop(gop), 
  mGep(gep) {
}


ImplFragmentorIterative::~ImplFragmentorIterative() 
{
  debug_func_cerr(5);

}

ImplFragmentorIterative::ImplFragmentorIterative( const ImplFragmentorIterative & src ) : 
    ImplFragmentor(src), 
    mDots(src.mDots), 
    mMinScore( src.mMinScore ),
    mGop(src.mGop),
    mGep(src.mGep) {
}

//------------------------------------------------------------------------------------------------
void ImplFragmentorIterative::performFragmentation( const Alignandum * row, 
						    const Alignandum * col, 
						    const Alignata * sample) {
    
    Alignata * original_dots = mDots;
    // true, if mDots contains a copy of orignal alignata object. 
    // Make sure, you do not delete it.
    bool is_copy = false;

    while ( 1 ) {
	Alignator * dottor    = makeAlignatorPublishAlignata( mDots );
	Alignator * alignator = makeAlignatorDotsSquared( mGop, mGep, dottor );

	Alignata * result = sample->getNew();
	alignator->align( row, col, result);

#ifdef DEBUG
	cout << "starting alignment" << *mDots << endl;	
	cout << "result" << *result << endl;
#endif

	delete alignator;
	delete dottor;
	
	if (result->getScore() >= mMinScore) {
	  mFragments->push_back( result );

	  // delete dots from dot-plot. Delete all dots in region      
	  Alignata * copy = makeAlignataMatrixUnsorted();
	  
	  copyAlignataRemoveRegion( copy, 
				    mDots, 
				    result->getRowFrom(),
				    result->getRowTo(),
				    result->getColFrom(),
				    result->getColTo());
	  
	  // delete old copy
	  if (is_copy) 
	    delete mDots;
	    
	  // substitute new copy instead of old copy
	  mDots = copy;
	  is_copy = true;
	    
	} else {
	    delete result;
	    break;
	}
	
    }

    if (is_copy) {
	delete mDots;
	mDots = original_dots;
    }

}
  

} // namespace alignlib




