//--------------------------------------------------------------------------------
// Project LibAlign
//
// Copyright (C) 2000 Andreas Heger All rights reserved
//
// Author: Andreas Heger <heger@ebi.ac.uk>
//
// $Id: Distor.h,v 1.1.1.1 2002/07/08 21:20:17 heger Exp $
//--------------------------------------------------------------------------------


#if HAVE_CONFIG_H
#include <config.h>
#endif

#ifndef DISTOR_H
#define DISTOR_H 1

#include "alignlib.h"

#include <iosfwd>
#include <string>

/**
   base class for methods calculating distance matrices from
   multiple alignments of protein sequences. If you have a set
   of single pairswise alignments, do it yourself :=).

   About the matrix:
   The matrix should be symmetric and support two commands:
	1. nrows():		for retrieving the number of rows
	2. operator(row,col):	for setting/retrieving values

   @author Andreas Heger
   @version $Id: Distor.h,v 1.1.1.1 2002/07/08 21:20:17 heger Exp $
   @short base class for calculating distance matrices from sequences
 */ 


namespace alignlib 
{
class MultipleAlignment;

class PhyloMatrix;

class Distor 
{
	// class member functions
	friend std::ostream & operator<<( std::ostream &, const Distor &);
	friend std::istream & operator>>( std::istream &, Distor &);
	
public:
	// constructors and desctructors
	Distor  ();

	Distor  (const Distor &);

	virtual ~Distor ();

	/** calculate a distance matrix from protein sequences
	@param multali multiple alignment of protein sequences
	@param matrix  matrix to use. If not supplied, the most basic matrix type will be used.
	 */
	virtual PhyloMatrix * calculateMatrix( PhyloMatrix * dest, const alignlib::MultipleAlignment * mali ) const = 0;

	/** return the maximum possible distance than can be achieved between two sequences */
	virtual PhyloMatrixValue getMaximumPossibleDistance() const = 0;

	/** Calculate distance between two rows from multiple alignment */
	virtual PhyloMatrixValue calculateDistance( const std::string & s_row_1, const std::string & s_row_2) const = 0;

};


}

#endif /* DISTOR_H */

