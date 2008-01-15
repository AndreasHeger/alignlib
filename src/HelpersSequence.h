/*
  alignlib - a library for aligning protein sequences

  $Id: HelpersSequence.h,v 1.2 2004/01/07 14:35:33 aheger Exp $

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

#ifndef HELPERS_SEQUENCE_H
#define HELPERS_SEQUENCE_H 1

#include <string>

#include "alignlib_fwd.h"
#include "alignlib_fwd.h"

namespace alignlib 
{
    
    /** Helper functions for class Alignment:
	
	1. factory functions
	
	2. accessor functions for default objects
	
	3. convenience functions
    */

    /* -----------------------------------------------------------------------------------------*/
    /* 1. factory functions */
    /** create a sequence from a NULL-terminated string */
	HAlignandum makeSequence( const char * sequence ); 
	
    HAlignandum makeSequence( const char * sequence, 
    		const HTranslator & translator );

    /** create a sequence from a string */
    HAlignandum makeSequence( const std::string & sequence );
   
    HAlignandum makeSequence( const std::string & sequence,
    		const HTranslator & translator );

    /** mutate a sequence according to a substitution matrix 
     * 
     * Initializes random generator with seed, if seed > 0.
     * */
    HAlignandum makeMutatedSequence( 
    			HAlignandum src, 
    			const HMutationMatrix & matrix,
    			const long seed = 0);
    
}

#endif	/* HELPERS_SEQUENCE_H */
