/*
  alignlib - a library for aligning protein sequences

  $Id: HelpersAlignment.h,v 1.3 2004/10/14 23:34:09 aheger Exp $

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

#ifndef HELPERS_ALIGNATA_H
#define HELPERS_ALIGNATA_H 1

#include <iosfwd>
#include <string>
#include "alignlib_fwd.h"
#include "alignlib_fwd.h"

namespace alignlib 
{

/** Helper functions for class Alignment:

	1. factory functions

	2. accessor functions for default objects

	3. convenience functions
 */

/* -------------------------------------------------------------------------------------------------------------------- */
/* 1. factory functions */
/** return a pointer to an Alignment-object */
HAlignment makeAlignmentSet(); 

/** return a pointer to an Alignment-object (iterate column-wise)*/
HAlignment makeAlignmentSetCol(); 

/** return a pointer to an Alignment-object */
HAlignment makeAlignmentHash(); 

/** return a pointer to an Alignment-object */
HAlignment makeAlignmentHashDiagonal(); 

/** return a pointer to an Alignment-object */
HAlignment makeAlignmentVector(); 

/** return a pointer to a new Alignment-object */
HAlignment makeAlignmentMatrixRow( long ndots = 0);

/** return a pointer to a new Alignment-object */
HAlignment makeAlignmentMatrixUnsorted( long ndots = 0);

/** return a pointer to a new Alignment-object */
HAlignment makeAlignmentMatrixDiagonal( long ndots = 0);


/* -------------------------------------------------------------------------------------------------------------------- */
/* 2. accessor functions for default objects */


/* -------------------------------------------------------------------------------------------------------------------- */
/* 3. convenience functions */

/** check if two alignmetns are identical
 * 
 * This method iterates over both alignments and checks if the
 * the same coordinates are returned. Thus, alignments that are 
 * sorted differently are not identical.
 */
bool checkAlignmentIdentity( 
		const HAlignment &,
		const HAlignment &,
		const bool invert = false );

/** enum describing the ways that two alignments can be combined
 * R: row
 * C: col
 */
typedef enum { RR, RC, CR, CC } CombinationMode;

/** return a pointer to an Alignment-object, that has been created by combining two others */
void combineAlignment( 
		HAlignment & dest, 
		const HAlignment & src1, 
		const HAlignment & src2, 
		const CombinationMode mode);

/** copy one alignment into another. I can not do this using a copy constructor, since virtual functions are not
     resolved by then. I also do not want a clone, for example if I want to change the underlying implementation

     @param row_from	beginning of segment to use
     @param row_to	end of segment to use
     @param col_from	beginning of segment to use
     @param col_to	end of segment to use
     @param diagonal_form beginning of tube to use
     @param diagonal_to end of tube to use
 */
void copyAlignment( HAlignment & dest, 
		const HAlignment & src, 
		Position row_from = NO_POS,
		Position row_to = NO_POS,
		Position col_from = NO_POS,
		Position col_to = NO_POS,
		Diagonal diagonal_from = -MAX_DIAGONAL,
		Diagonal diagonal_to = MAX_DIAGONAL
);

/** copy one alignment into another. In contrast to copyAlignment, this function keeps the dots outside
     of region.
     @param row_from	beginning of segment to use
     @param row_to	end of segment to use
     @param col_from	beginning of segment to use
     @param col_to	end of segment to use
     @param diagonal_form beginning of tube to use
     @param diagonal_to end of tube to use
 */
void copyAlignmentRemoveRegion( HAlignment & dest, 
		const HAlignment & src, 
		Position row_from = NO_POS,
		Position row_to = NO_POS,
		Position col_from = NO_POS,
		Position col_to = NO_POS,
		Diagonal diagonal_from = 1,
		Diagonal diagonal_to = 0
);

/** copy one alignment into another using filter 
 */
void copyAlignment( HAlignment & dest, 
		const HAlignment & src, 
		const HAlignment & filter, 
		const CombinationMode mode);

/** add one alignment to another
 */
void addAlignment2Alignment( HAlignment & dest, const HAlignment & src );

/** add one alignment to another. Map src using map_src2new.
 */
void addMappedAlignment2Alignment( HAlignment & dest, 
		const HAlignment & src, 
		const HAlignment & map_src2new,
		const CombinationMode mode );

/** add one alignment to another. Map both row and column.
 */
void addMappedAlignments2Alignment( HAlignment & dest, 
		const HAlignment & src, 
		const HAlignment & map_src_row2dest_row, 
		const HAlignment & map_src_col2dest_col );


/** create an identity alignment between residues from and to in row using an offset for col */
void fillAlignmentIdentity( 
		HAlignment & dest, 
		Position row_from, 
		Position row_to, 
		Position col_offset = 0);

/** fill gaps in an alignment by doing a local alignment in each
     region.
 */
void fillAlignmentGaps( 
		HAlignment & dest,
		const HAlignator & alignator,
		const HAlignandum & row,
		const HAlignandum & col );

/** remove residues from an alignment, that are part of another alignment
     @param mode: specifies, which residues are looked up. If mode = RR, then every pairs is eliminated from dest,
     where the row is also present as a row-residue in filter.
 */
void filterAlignmentRemovePairs( HAlignment & dest, 
		const HAlignment & filter, 
		const CombinationMode mode );

/** remove residues from an alignment, that are part of another alignment
     @param mode: specifies, which residues are looked up. If mode = RR, then every pairs is eliminated from dest,
     where the row is also present as a row-residue in filter.
 */
void filterAlignmentRemovePairwiseSorted( HAlignment & dest, 
		const HAlignment & filter, 
		const CombinationMode mode );

/** rescore residues in an alignment. 
 */
void rescoreAlignment( 
		HAlignment & dest,
		const HAlignandum & row,
		const HAlignandum & col,
		const HScorer & scorer );

/** rescore residues in an alignment setting them all to the same score. 
 */
void rescoreAlignment( 
		HAlignment & dest,
		const Score score = 0);

/** calculate Alignment score given gap-penalties for row and column */
void calculateAffineScore( HAlignment & dest, 
		const Score gop, 
		const Score gep );


/** fill an alignment with a repeat unit from a wrap-around alignment */
void fillAlignmentRepeatUnit( 
		HAlignment & dest, 
		const HAlignment & source,
		const Position first_row_residue = NO_POS,
		const bool skip_negative_ends = false);

/** 
    return the maps of row/col of an alignment to the summation of the alignment. This is useful for 
    building multiple alignemnts.

    @param	map_row2combined     
    @param	map_col2combined
    @param      source
    @param	insert_gaps_row
    @param	insert_gaps_col
    @param	use_end_row
    @param	use_end_col
    @param	row_length
    @param	col_length

 */

void fillAlignmentSummation( HAlignment & dest1, 
		HAlignment & dest2, 
		const HAlignment & src,
		const bool insert_gaps_row = true,
		const bool insert_gaps_col = true, 
		const bool use_end_row = false,
		const bool use_end_col = false, 
		const Position row_length = NO_POS,
		const Position col_length = NO_POS);


/** complement a pairwise alignment. If there is a gap of the same length in both row and
     col, the corresponding residues are added to the alignment.
 */
void complementAlignment( 
		HAlignment & dest, 
		const Position max_length );

/** remove all those residues from an alignmnent, which are not
     in sequential. 
     
     This ensures, that col_i < col_i+1 and row < row_i+1
     
     Useful with AlignmentVector.
 */
void flattenAlignment( HAlignment & dest );

/** split an alignment, if there are gaps larger than a certain 
 * 	threshold either in row or col or both.
 */
HFragmentVector splitAlignment( 
		const HAlignment & src, 
		const int max_gap_width,
		bool split_row = true,
		bool split_col = true);

/** split an alignment at points of intersection with another alignment.
 */ 
HFragmentVector splitAlignment( 
		const HAlignment & src1, 
		const HAlignment & src2, 
		const CombinationMode mode );

/** starting from the ends of an alignment, remove 
    residues which do not contribute to a positive score.
 */
void pruneAlignment( 
		HAlignment & src,
		const Score gop,
		const Score gep);

/** calculate percent similarity of alignment */
double calculatePercentSimilarity( 
		const HAlignment & src);  

/** calculate percent identity of alignment. 
 */
double calculatePercentIdentity (
		const HAlignment & src, 
		const HAlignandum & row, 
		const HAlignandum & col); 

/** remove small fragments from alignment.
 * 
 * This method removes fragments from an alignment. A fragment
 * is a part of an alignment that is short (max_fragment_length)
 * and surrounded by large gaps (min_gap_length).
 */
void removeFragments( 
		HAlignment & dest,
		const unsigned int window_length,
		const unsigned int min_gap_length,
		const Position row_length = NO_POS);

}

#endif	/* HELPERS_ALIGNATA_H */






