/*
  alignlib - a library for aligning protein sequences

  $Id: HelpersMultipleAlignment.h,v 1.5 2004/03/19 18:23:40 aheger Exp $

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

#ifndef HELPERS_MULTIPLE_ALIGNMENT_H
#define HELPERS_MULTIPLE_ALIGNMENT_H 1

#include <vector>
#include "alignlib_fwd.h"
#include "MultipleAlignment.h"

namespace alignlib 
{

    template<class T> class Matrix;
    typedef Matrix<unsigned int> CountsMatrix;
    typedef std::vector<double> VectorDouble;

    /** Helper functions for class Alignment:
	
	1. factory functions
	
	2. accessor functions for default objects
	
	3. convenience functions
    */
    
    /* -------------------------------------------------------------------------------------------------------------------- */
    /* 1. factory functions */

    /** return empty alignment */
    HMultipleAlignment makeMultipleAlignment(); 

    /** return empty alignment */
    HMultipleAlignment makeMultipleAlignmentDots( bool compress_unaligend_columns = true,
						   int max_insertion_length = -1);

    /* -------------------------------------------------------------------------------------------------------------------- */
    /* 2. accessor functions for default objects */
    

    /* -------------------------------------------------------------------------------------------------------------------- */
    /* 3. convenience functions */

    /** fill multiple alignment ali from contents of a file */
    /*
    HMultipleAlignment fillMultipleAlignment( HMultipleAlignment ali,
					       const char * filename);
	*/
    /** fill multiple alignment ali from memory, sequences contains the concatenated rows without separator,
	however there has to be a \0 at the end of sequences
	
	@param ali		alignment to store data in.
	@param sequences	concatenated sequences to use
	@param nsequences	number of sequences.
    */
    void fillMultipleAlignment( 
    		HMultipleAlignment & ali,
    		const std::string & sequences,
    		int nsequences);
    
    /** fill multiple alignment ali from contents of a file, use an Alignatum-object for parsing
     */
    /* HMultipleAlignment fillMultipleAlignment( HMultipleAlignment ali, 
					       const char * filename, 
					       const Alignatum * alignatum_template);
					       */
    /** fill multiple alignment ali from contents of a file, use an Alignatum-object for parsing */
    /*
    HMultipleAlignment fillMultipleAlignment( HMultipleAlignment ali, 
					       const char * filename, 
					       const Alignatum * alignatum_template);
*/
    /** extract a multiple Alignment object from a stream in FASTA format */
    /*
    HMultipleAlignment extractMultipleAlignmentFasta( HMultipleAlignment ali, 
						       std::istream & input );
    */
    /** get Convservation-string for multiple alignment. This returns a string, where 
	each residue is marked, which is conserved at least > cutoff % */        
    /*
    std::string calculateConservation( const HMultipleAlignment & mali, 
    		Frequency min_frequency);
	*/

    /** calculate counts in mali categorised. The first row in Matrix is empty, so that
	the numbering is consistent with the numbering of columns in the multiple alignment.
	Only include num_rows rows in the alignment. If num_rows is not set (0), all
	rows are taken.
    */
    /*
    CountsMatrix * makeCountsByCategory( const HMultipleAlignment & mali, 
					 const unsigned int * map_residue2category = NULL);
	*/
    /** make a map from residues to categories. The following order has been suggested by Hannes for
	surface area calculations:
          'K': 1, 'R': 1,
	'D': 2, 'E': 2,
	'H': 3, 'F': 3, 'W':3, 'Y': 3, 'C': 3,
	'N': 4, 'Q': 4, 'S':4, 'T': 4,
	'A': 5, 'I': 5, 'L': 5,'M': 5, 'V': 5,
	'G': 6, 'P': 6
    */
    const unsigned int * getMapResidue2CategorySurface();

    /** mapping, where each residue is mapped to its own 
	category 
    */
    const unsigned int * getMapResidue2CategoryAll();


    /** return a vector of entropies calculated for a CountsMatrix
     */
    VectorDouble * makeEntropyVector( const CountsMatrix * src);

    /** copy a multiple alignment
     */
    void copyMultipleAlignment( 
    		HMultipleAlignment & dest, 
    		const HMultipleAlignment & src,
    		unsigned int first_row = 0,
    		unsigned int last_row = 0 );
}

#endif	/* HELPERS_MULTIPLE_ALIGNMENT_H */
