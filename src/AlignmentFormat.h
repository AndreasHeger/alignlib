#ifndef ALIGNMENTFORMAT_H_
#define ALIGNMENTFORMAT_H_

//--------------------------------------------------------------------------------
// Project alignlib
//
// Copyright (C) 2000 Andreas Heger All rights reserved
//
// Author: Andreas Heger <heger@ebi.ac.uk>
//
// $Id$
//--------------------------------------------------------------------------------

/** 
 * Alignment data structures. These classes convert various alignment formats
 * to/from alignment objects. Member data is publich for easy access to your
 * own formatting functions.
 * 
 * Usage example:
 * 
 * HAligment ali = makeAlignmentVector();
 * ali->fill(...);
 * 
 * // write blocks
 * AlignmentFormatBlocks blocks_out(ali);
 * 
 * std::cout << ali << std::endl; 
 * // read blocks
 * AlignmentFormatBlocks blocks_in;
 * std::cin >> blocks_in >> std::endl;
 * blocks_in.copy( ali );
 * 
 */

#if HAVE_CONFIG_H
#include <config.h>
#endif

#include <iosfwd>
#include <string>

#include "alignlib_fwd.h"

namespace alignlib 
{

/** base class - keeps track of alignment coordinates
 */
struct AlignmentFormat
{
	
	// constructors and desctructors
	AlignmentFormat ();

	AlignmentFormat( const AlignmentFormat &);

	virtual ~AlignmentFormat();

	/** fill format from alignment
		@param src Alignment to parse
	 */
	virtual void fill( const HAlignment & src);

	/** fill alignment from format
	 * 	@param dest Alignment 
	 */
	virtual void copy( HAlignment & dest ) const;
	
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
};

/**
	Data structures for aligned blocks alignment format
	
	This is a convenience structure for importing/exporting
	pairwise alignments.
	
   	@author Andreas Heger
   	@version $Id$
   	@short Data structure of aligned blocks.
 
*/ 
struct AlignmentFormatBlocks : public AlignmentFormat
{
	// class member functions
	friend std::ostream & operator<<( std::ostream &, const AlignmentFormatBlocks &);
	friend std::istream & operator>>( std::istream &, AlignmentFormatBlocks &);
	
	// constructors and desctructors
	AlignmentFormatBlocks ();

	AlignmentFormatBlocks( const HAlignment & src);
	
	AlignmentFormatBlocks (const AlignmentFormatBlocks &);

	virtual ~AlignmentFormatBlocks ();

	/** fill blocks from alignment
		@param src Alignment to parse
	 */
	virtual void fill( const HAlignment & src);

	/** fill Alignment object with blocks
	 * 	@param dest Alignment 
	 */
	virtual void copy( HAlignment & dest ) const;
	
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

/**
	Data structures for "Emissions" alignment format
	
	This is a convenience structure for importing/exporting
	pairwise alignments.
	
	This format stores the alignment in two strings for 
	row and col, respectively. Each string represents
	the emissions (prefixed by +) and insertions 
	(prefixed by -) for the alignment.
	
   	@author Andreas Heger
   	@version $Id$
   	@short Data structure of emissions alignment format.
 
*/ 

struct AlignmentFormatEmissions : public AlignmentFormat
{
	// class member functions
	friend std::ostream & operator<<( std::ostream &, const AlignmentFormatEmissions &);
	friend std::istream & operator>>( std::istream &, AlignmentFormatEmissions &);
	
	// constructors and desctructors
	AlignmentFormatEmissions ();

	AlignmentFormatEmissions( const HAlignment & src);
	
	AlignmentFormatEmissions (const AlignmentFormatEmissions &);

	virtual ~AlignmentFormatEmissions ();

	/** fill blocks from alignment
		@param src Alignment to parse
	 */
	virtual void fill( const HAlignment & src);

	/** fill Alignment object with blocks
	 * 	@param dest Alignment 
	 */
	virtual void copy( HAlignment & dest ) const;
	
	/** the alignment for row 
	 * 
	 */
	std::string mRowAlignment;

	/** the alignment for col
	 * 
	 */
	std::string mColAlignment;
	
};

/**
	Data structures for "Explicit" alignment format
	
	This is a convenience structure for importing/exporting
	pairwise alignments.

	This format represents the alignment as two strings
	for row and col, respectively. The alignment is 
	represented as aligned characters.
	
   	@author Andreas Heger
   	@version $Id$
   	@short Data structure of explicit alignment format.
 
*/ 

struct AlignmentFormatExplicit : public AlignmentFormat
{
	// class member functions
	friend std::ostream & operator<<( std::ostream &, const AlignmentFormatExplicit &);
	friend std::istream & operator>>( std::istream &, AlignmentFormatExplicit &);
	
	// constructors and desctructors
	AlignmentFormatExplicit ();

	AlignmentFormatExplicit( const HAlignment & src,
			const HAlignandum & row,
			const HAlignandum & col);
	
	AlignmentFormatExplicit (const AlignmentFormatExplicit &);

	virtual ~AlignmentFormatExplicit ();

	/** fill blocks from alignment
		@param src Alignment to parse
	 */
	virtual void fill( const HAlignment & src,
			const HAlignandum & row,
			const HAlignandum & col);

	/** fill Alignment object with blocks
	 * 	@param dest Alignment 
	 */
	virtual void copy( HAlignment & dest) const;
	
	/** the alignment for row 
	 * 
	 */
	std::string mRowAlignment;

	/** the alignment for col
	 * 
	 */
	std::string mColAlignment;
	
};

/**
	Data structures for "Diagonals" alignment format
	
	This is a convenience structure for importing/exporting
	pairwise alignments.
	
	The diagonals alignment stores the alignment a single
	string. The string records emissions and insertions
	along each digaonal. This alignment format is useful for
	storing dotplots.

    Although any alignment class can be written in this format, 
    it is best to use for those that are sorted by diagonal, for
    example, MatrixDiagonal. This achieves the highest compression.
     	
   	@author Andreas Heger
   	@version $Id$
   	@short Data structure of "Diagonals" alignment format.
 
*/ 

struct AlignmentFormatDiagonals : public AlignmentFormat
{
	// class member functions
	friend std::ostream & operator<<( std::ostream &, const AlignmentFormatDiagonals &);
	friend std::istream & operator>>( std::istream &, AlignmentFormatDiagonals &);
	
	// constructors and desctructors
	AlignmentFormatDiagonals ();

	AlignmentFormatDiagonals( const HAlignment & src);
	
	AlignmentFormatDiagonals (const AlignmentFormatDiagonals &);

	virtual ~AlignmentFormatDiagonals ();

	/** fill format from alignment

	 	 @param src Alignment to parse
	 	 @param reverse		reverse column and row 
	 	 @param row_from	beginning of segment to use
	     @param row_to	end of segment to use
	     @param col_from	beginning of segment to use
	     @param col_to	end of segment to use
	     @param diagonal_form beginning of tube to use
	     @param diagonal_to end of tube to use

		 The alignment can be restricted to a region specifying the columns. 
		 A further filter can be applied, that only saves a band. 

	     If diagonal_from is larger than diagonal_to, then the whole range is used.
	     
    The diagonal format is:

    diagonal:-n1+n2-n3+n4;diagonal;...
	     
	 */
	virtual void fill( 
			const HAlignment & src,
			const bool reverse,
			const Position row_from = NO_POS,
			const Position row_to = NO_POS,
			const Position col_from = NO_POS,
			const Position col_to = NO_POS,
			const Diagonal diagonal_from = MAX_DIAGONAL,
			const Diagonal diagonal_to = -MAX_DIAGONAL
	);

	virtual void fill( 
			const HAlignment & src );
	
	/** fill Alignment object from format
		@param dest 	Alignment to fill
		@param reverse 	reverse column and row. 
	 */
	virtual void copy(
			HAlignment & dest,
			const bool reverse) const;	

	virtual void copy( HAlignment & dest ) const;
		
	/** the alignment 
	 * 
	 */
	std::string mAlignment;
	
};



}

#endif
