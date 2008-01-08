/*
  alignlib - a library for aligning protein sequences

  $Id: HelpersIterator2D.h,v 1.2 2004/01/07 14:35:32 aheger Exp $

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

#ifndef HELPERS_SCORER_H
#define HELPERS_SCORER_H 1

#include <iosfwd>

#include "alignlib.h"
#include "alignlib_fwd.h"
#include "alignlib_default.h"
#include "Scorer.h"

namespace alignlib 
{
    
    /** Helper functions for class Alignment:
	
	1. factory functions
	
	2. accessor functions for default objects
	
	3. convenience functions
    */


  /* 1. factory functions */
	HScorer makeScorer( const HAlignandum & row,
			const HAlignandum & col );
  
	/** return a placeholder scorer */
	HScorer makeScorer();
	
  /* 2. accessor functions for default objects */

  /* 3. convenience functions */
  
  DEFINE_DEFAULT( HScorer, getDefaultScorer, setDefaultScorer );
}

#endif	/* HELPERS_SCORER_H */
