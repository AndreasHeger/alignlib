/*
  alignlib - a library for aligning protein sequences

  $Id: HelpersLogOddor.h,v 1.2 2004/01/07 14:35:32 aheger Exp $

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

#ifndef HELPERS_LOGODDOR_H
#define HELPERS_LOGODDOR_H 1

#include "alignlib.h"
#include "LogOddor.h"

namespace alignlib {
    
    /** Helper functions for class Alignata:
	
	1. factory functions
	
	2. accessor functions for default objects
	
	3. convenience functions
    */

    
    /* -------------------------------------------------------------------------------------------------------------------- */
    /* 1. factory functions */
    
    const LogOddor * makeLogOddorUniform( Score scale_factor = 1.0);

    const LogOddor * makeLogOddorDirichlet( Score scale_factor = 1.0);

    /* -------------------------------------------------------------------------------------------------------------------- */
    /* 2. accessor functions for default objects */

    /** gets the default LogOddor object */ 
    const LogOddor * getDefaultLogOddor();

    /** sets the default LogOddor object */
    const LogOddor * setDefaultLogOddor( const LogOddor * logoddor);
    

    /* -------------------------------------------------------------------------------------------------------------------- */
    /* 3. convenience functions */
    

}

#endif	/* HELPERS_LOGODDOR_H */
