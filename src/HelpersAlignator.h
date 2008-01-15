/*
  alignlib - a library for aligning protein sequences

  $Id: HelpersAlignator.h,v 1.4 2005/02/24 11:07:25 aheger Exp $

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

#ifndef HELPERS_ALIGNATOR_H
#define HELPERS_ALIGNATOR_H 1

#include "alignlib_fwd.h"
#include "alignlib_fwd.h"

namespace alignlib 
{
    
    /** Helper functions for class Alignment:
	
	1. factory functions
	
	2. accessor functions for default objects
	
	3. convenience functions
    */

  /* 1. factory functions */

    /** Perform full dynamic programming  */
	HAlignator makeAlignatorDPFull( 
			AlignmentType alignment_type,
			Score gop, Score gep, 
			bool penalize_left = false, 
			bool penalize_right = false );

    /** make an alignator object, which aligns identical residues */
    HAlignator makeAlignatorIdentity();

    /** make an alignator object, which aligns identical residues */
    HAlignator makeAlignatorSimilarity();

    /** make an alignator object, which aligns similar tuples */
    HAlignator makeAlignatorTuples(int ktuple = 3 );

    /** make an alignator object, which returns a dummy alignments */
    HAlignator makeAlignatorDummy( const HAlignment & ali );

    /** make an alignator object, which returns a dummy alignment */
    HAlignator makeAlignatorPublishAlignment( HAlignment & ali );
    
    /** make an alignator object, which does a dot-alignment with wrapping around. */
    HAlignator makeAlignatorDotsWrap(
    		const HAlignator & alignator,    		
    		Score gop, 
    		Score gep );

    /** make an alignator object, which does a dot-alignment with wrapping around. */
    HAlignator makeAlignatorDotsSquared(
    		const HAlignator & alignator,
    		Score gop, 
    		Score gep );
    
    /** make an alignator object, which does a dot-alignment with wrapping around. */
    HAlignator makeAlignatorDotsSquaredDiagonal(
    		const HAlignator & alignator,     		
    		Score gop, 
    		Score gep, 
    		Score diagnal_gop = 0,
    		Score diagonal_gep = 0 );

    /** make an alignator object, which aligns fragments. */
    HAlignator makeAlignatorFragmentsSquared(
    		Score gop, 
    		Score gep, 
    		const HFragmentor & fragmentor );

    /** alignator object for iterative alignment
     * 
     * Aligns two Alignandum objects iteratively using a template
     * alignator object until the alignment score drops below @min_score. 
     * The template alignator object is copied.
     */
    HAlignator makeAlignatorIterative( 
    		const HAlignator & alignator, 
    		Score min_score);
    
    // compatibility functions
    inline HAlignator makeFullDP( 
    		Score gop, Score gep, 
    		bool penalize_left = false, 
    		bool penalize_right = false )
      {
    	return makeAlignatorDPFull( ALIGNMENT_LOCAL,
    			gop, gep, false, false);
      }
    

    /* 2. accessor functions for default objects */
    

    /* 3. convenience functions */
        
}

#endif	/* HELPERS_ALIGNATOR_H */
