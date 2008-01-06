/*
  alignlib - a library for aligning protein sequences

  $Id: HelpersRegularizor.h,v 1.2 2004/01/07 14:35:32 aheger Exp $

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

#ifndef HELPERS_REGULARIZOR_H
#define HELPERS_REGULARIZOR_H 1

#include "alignlib.h"

namespace alignlib 
{

	class Regularizor;

    /** Helper functions for class Alignment:
	
	1. factory functions
	
	2. accessor functions for default objects
	
	3. convenience functions
    */

    
    /* -------------------------------------------------------------------------------------------------------------------- */
    /* 1. factory functions */

	/** empty regularizor - computes frequencies from counts */
    Regularizor * makeRegularizor();
    
    /** regularizor according to Tatusov et al.()
     * for own parameterization. 
     */
    Regularizor * makeRegularizorTatusov( const SubstitutionMatrix * matrix,
    		const FrequencyVector & background,
    		const double & beta, 
    		const double & lambda ); 

    /** regularizor according to Tatusov et al.() parameterized according to PSIBLAST
     */    
    Regularizor * makeRegularizorPsiblast();
    
    Regularizor * makeRegularizorDirichlet( Count fade_cutoff = 0);
    
    Regularizor * makeRegularizorDirichletHash( Count fade_cutoff = 0);
    
    Regularizor * makeRegularizorDirichletInterpolate( Count fade_cutoff = 0);
    
    Regularizor * makeRegularizorDirichletPrecomputed( Count fade_cutoff = 0);

    /* -------------------------------------------------------------------------------------------------------------------- */
    /* 2. accessor functions for default objects */

    /** gets the default Regularizor object */ 
    const Regularizor * getDefaultRegularizor();

    /** sets the default Regularizor object 
     * The library obtains ownership of the regularizor
     * */
    void setDefaultRegularizor( Regularizor * regularizor);

    /* -------------------------------------------------------------------------------------------------------------------- */
    /* 3. convenience functions */
    

}

#endif	/* HELPERS_REGULARIZOR_H */
