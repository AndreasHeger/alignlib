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

#include "alignlib.h"

namespace alignlib {
    
    /** Helper functions for class Alignment:
	
	1. factory functions
	
	2. accessor functions for default objects
	
	3. convenience functions
    */

    
    class Alignandum;
    class SubstitutionMatrix;

    /* defintion of matrix for mutating a sequence */
    template<class T> class Matrix;
    typedef Matrix<double> MutationMatrix;

    /* -----------------------------------------------------------------------------------------*/
    /* 1. factory functions */
    /** create a sequence from a NULL-terminated string */
    Alignandum * makeSequence( const char * sequence );

    /** create a sequence from a string */
    Alignandum * makeSequence( const std::string & sequence );

    /** mutate a sequence according to a substitution matrix */
    Alignandum * makeMutatedSequence( Alignandum * src, 
				      const MutationMatrix * matrix );

    /* ----------------------------------------------------------------------------------------------*/
    /* 2. accessor functions for default objects */
    

    /* ----------------------------------------------------------------------------------------------*/ 
    /* 3. convenience functions */
    /** create a sequence from a stream */
    Alignandum * extractSequence( std::istream & input );

    /** extract a sequence in Fasta-Format from a stream */
    Alignandum * extractSequenceFasta( std::istream & input, std::string & description );

    /** create a sequence from a filename */
    Alignandum * readSequence( const char * filename );

    /** set random seed */
    void SetRandomSeed( const long seed );
    
}

#endif	/* HELPERS_SEQUENCE_H */
