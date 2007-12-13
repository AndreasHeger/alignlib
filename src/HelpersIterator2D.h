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

#ifndef HELPERS_ITERATOR2D_H
#define HELPERS_ITERATOR2D_H 1

#include <iosfwd>

#include "alignlib.h"
#include "Iterator2D.h"

namespace alignlib {
    
    /** Helper functions for class Alignment:
	
	1. factory functions
	
	2. accessor functions for default objects
	
	3. convenience functions
    */


  class Alignandum;
 
  /* 1. factory functions */
  Iterator2D * makeIterator2DFull( const Alignandum * row = NULL,
				   const Alignandum * col = NULL );
  
  Iterator2D * makeIterator2DBanded( const Alignandum * row = NULL,
				     const Alignandum * col = NULL,
				     const Diagonal lower_diagonal = 0,
				     const Diagonal upper_diagonal = 0);

  /* 2. accessor functions for default objects */
  const Iterator2D * getDefaultIterator2D();
  
  /** set the default iterator.
   * 
   * The library takes ownership of the supplied iterator
   */
  void setDefaultIterator2D( Iterator2D *);

  /* 3. convenience functions */
}

#endif	/* HELPERS_ITERATOR2D_H */
