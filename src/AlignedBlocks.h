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

#ifndef ALIGNED_BLOCKS_H
#define ALIGNED_BOCKS_H 1

#include <iosfwd>
#include <string>

#include "alignlib.h"

/**
	Data structure for aligned blocks.
	
	This is a convenience structure for importing/exporting
	pairwise alignments.
	
   	@author Andreas Heger
   	@version $Id$
   	@short Data structure of aligned blocks.
 
*/ 


namespace alignlib 
{
struct AlignedBlocks
{
	// class member functions
	friend std::ostream & operator<<( std::ostream &, const AlignedBlocks &);
	friend std::istream & operator>>( std::istream &, AlignedBlocks &);
	
	// constructors and desctructors
	AlignedBlocks  ( const Alignata * src = NULL);

	AlignedBlocks  (const AlignedBlocks &);

	virtual ~AlignedBlocks ();

	/** fill blocks from alignment
	@param src Alignment to parse
	 */
	void fill( const Alignata * src);

	/** fill Alignata object with blocks
	 * @param dest Alignment 
	 */
	Alignata * copy( Alignata * dest ) const;
	
	/** start of aligned blocks in row*/
	Position mRowFrom;
	
	/** end of aligned blocks
	 */
	Position mRowTo;
	
	/** start of aligned blocks in col 
	 */
	Position mColFrom;
	
	/** end of aligned blocks in col  
	 */
	Position mColTo;
	
	/** vector with starts of aligned blocks. 
	 * 
	 * The start positions are relative to the offset in @mRowFrom */
	PositionVector mRowStarts;

	/** vector with ends of aligned blocks.
	 * 
	 * The start positions are relative to the offset in @mColFrom */
	PositionVector mColStarts;
	
	/** vector with block sizes */
	PositionVector mBlockSizes;
	
	
};


}

#endif /* AlignedBlocks_H */

