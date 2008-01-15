/*
  alignlib - a library for aligning protein sequences

  $Id: HelpersAlignatum.h,v 1.2 2004/01/07 14:35:32 aheger Exp $

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

#ifndef HELPERS_ALIGNATUM_H
#define HELPERS_ALIGNATUM_H 1

#include "alignlib_fwd.h"
#include "alignlib_fwd.h"

namespace alignlib 
{
    
    /** Helper functions for class Alignment:
	
	1. factory functions
	
	2. accessor functions for default objects
	
	3. convenience functions
    */

    /* -------------------------------------------------------------------------------------------------------------------- */
    /* 1. factory functions */
    HAlignatum makeAlignatum( 
    		const std::string & src, 
    		const Position from = NO_POS, 
    		const Position to = NO_POS);
    
    HAlignatum makeAlignatum(
    		const HAlignandum & src, 
    		const HAlignment & map_this2new,
    		const Position max_length = 0);

    HAlignatum makeAlignatum(
    		const HAlignandum & src,
    		const Position from = NO_POS,
    		const Position to = NO_POS);

    HAlignatum makeAlignatumFasta(
    		const std::string & description,
    		const std::string & src,
    		const Position from = 0, 
    		const Position to = NO_POS);
    
    /** return an empty Alignatum object, that 
	can store a description */
    HAlignatum makeAlignatumFasta();

    /** return an empty Alignatum object, that 
	can be filled for example by the function Read */
    HAlignatum makeAlignatum();

    /* -------------------------------------------------------------------------------------------------------------------- */
    /* 2. accessor functions for default objects */
    

    /* -------------------------------------------------------------------------------------------------------------------- */
    /* 3. convenience functions */
    

}

#endif	/* HELPERS_ALIGNATUM_H */
