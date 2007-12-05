/*
  alignlib - a library for aligning protein sequences

  $Id: HelpersSubstitutionMatrix.h,v 1.2 2004/01/07 14:35:33 aheger Exp $

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

#ifndef HELPERS_SUBSTITUTION_MATRIX_H
#define HELPERS_SUBSTITUTION_MATRIX_H 1

#include "alignlib.h"

namespace alignlib {
    
    /** Helper functions for class Alignata:
	
	1. factory functions
	
	2. accessor functions for default objects
	
	3. convenience functions
    */
    
    class SubstitutionMatrix;

    /* -------------------------------------------------------------------------------------------------------------------- */
    /* 1. factory functions */

    /** create a Substitution Matrix given a memory location of log-odds scores */
    SubstitutionMatrix * makeSubstitutionMatrixAA( ScoreColumn * matrix, bool this_own = false);

    /** create the identity substitution matrix. Identities score as 1, mismatches as -1.
     */
    SubstitutionMatrix * makeSubstitutionMatrixAAIdentity( const Score match = 10,
							   const Score mismatch = -1);

    /** read a simple substitution matrix adn return it. Do not forget to tell the substitution matrix
	to clear its contents afterwards
    */
    SubstitutionMatrix * readSubstitutionMatrixAA( const char * filename ); 

    /* -------------------------------------------------------------------------------------------------------------------- */
    /* 2. accessor functions for default objects */

    /** gets the default SubstitionMatrix.
     * 
     * The library keeps ownership of the object 
     * */ 
    const SubstitutionMatrix * getDefaultSubstitutionMatrix();
    
    /** sets the default SubstitutionMatrix.
     * 
     * The library takes ownership of the default matrix.
     * The return type is void because of constraints in py++/boost.python
     * */
    void setDefaultSubstitutionMatrix( SubstitutionMatrix * matrix);

    /* -------------------------------------------------------------------------------------------------------------------- */
    /* 3. convenience functions */

}

#endif	/* HELPERS_SUBSTITUTION_MATRIX_H */



