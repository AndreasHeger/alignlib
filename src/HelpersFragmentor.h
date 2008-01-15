/*
  alignlib - a library for aligning protein sequences

  $Id: HelpersFragmentor.h,v 1.2 2004/01/07 14:35:32 aheger Exp $

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


#if HAVE_CONFIG_H
#include <config.h>
#endif

#ifndef HELPERS_FRAGMENTOR_H
#define HELPERS_FRAGMENTOR_H 1

#include <iosfwd>

#include "alignlib_fwd.h"
#include "Fragmentor.h"

namespace alignlib {
    
    /** Helper functions for class Alignment:
	
	1. factory functions
	
	2. accessor functions for default objects
	
	3. convenience functions
    */


  class Alignator;
 
  /* 1. factory functions */
  HFragmentor makeFragmentorDiagonals(
		  const HAlignator & alignator,
		  Score gop, 
		  Score gep );
  
  HFragmentor makeFragmentorIterative( 
		  const HAlignment & dots, 
		  Score min_score,
		  Score gop,
		  Score gep );

  HFragmentor makeFragmentorRepetitive( 
		  const HAlignator & alignator,
		  Score min_score );

  /* 2. accessor functions for default objects */
    

  /* 3. convenience functions */
  void writeFragments( 
		  std::ostream & output,
		  const HFragmentVector & fragments);
 
  /** rescore fragments. 
      The score of a fragment is the number of aligned residues plus gap penalties. 
      Gaps are penalized using affine gap penalties.
  */
  void rescoreFragmentsNumberGaps( 
		  HFragmentVector & fragments, 
		  Score gop = 0, 
		  Score gep = 0);

}

#endif	/* HELPERS_FRAGMENTOR_H */
