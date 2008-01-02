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
#include "alignlib.h"
#include "SubstitutionMatrix.h"

namespace alignlib 
{

/** Helper functions for class Alignment:

	1. factory functions

	2. accessor functions for default objects

	3. convenience functions
 */

class Alignator;
class Alignment;
class Alignandum;

/* -------------------------------------------------------------------------------------------------------------------- */
/* 1. factory functions */
/** return a pointer to an Alignment-object */
Alignment * makeAlignmentSet(); 

/** return a pointer to an Alignment-object (iterate column-wise)*/
Alignment * makeAlignmentSetCol(); 

/** return a pointer to an Alignment-object */
Alignment * makeAlignmentHash(); 

/** return a pointer to an Alignment-object */
Alignment * makeAlignmentHashDiagonal(); 

/** return a pointer to an Alignment-object */
Alignment * makeAlignmentVector(); 

/** return a pointer to a new Alignment-object */
Alignment * makeAlignmentMatrixRow( long ndots = 0);

/** return a pointer to a new Alignment-object */
Alignment * makeAlignmentMatrixUnsorted( long ndots = 0);

/** return a pointer to a new Alignment-object */
Alignment * makeAlignmentMatrixDiagonal( long ndots = 0);


/* -------------------------------------------------------------------------------------------------------------------- */
/* 2. accessor functions for default objects */


/* -------------------------------------------------------------------------------------------------------------------- */
/* 3. convenience functions */


/** Print aligment in table format */
void writeAlignmentTable( 
		std::ostream & output, 
		const Alignment * src,
		unsigned int ncols = 8,
		bool with_scores = true);

/** write compressed alignment into stream. The alignment can be restricted to
     a region specifying the columns. If you want to restrict it to row-residues,
     you can get those boundaries by calling mapRowToCol

	 @param output 		output stream
     @param src	Alignment
     @param col_from	beginning of segment to use
     @param col_to		end of segment to use
 */
void writeAlignmentCompressed(
		std::ostream & output,
		const Alignment * src, 
		Position col_from = NO_POS,
		Position col_to = NO_POS
);

/** write compressed alignment into streams. The alignment can be restricted to
     a region specifying the columns. A further filter can be
     applied, that only saves a tube. 

     If diagonal_from is larger than diagonal_to, then the whole range is used.

	 @param output 		output stream
     @param src	Alignment
     @param reverse  whether to reverse row and column
     @param row_from	beginning of segment to use
     @param row_to	end of segment to use
     @param col_from	beginning of segment to use
     @param col_to	end of segment to use
     @param diagonal_form beginning of tube to use
     @param diagonal_to end of tube to use
 */
void writeAlignmentCompressedDiagonal(
		std::ostream & output,
		const Alignment * src, 
		bool reverse = false,
		Position row_from = NO_POS,
		Position row_to = NO_POS,
		Position col_from = NO_POS,
		Position col_to = NO_POS,
		Diagonal diagonal_from = MAX_DIAGONAL,
		Diagonal diagonal_to = -MAX_DIAGONAL
);

/** write an alignment in rsdb-format */
void writeAlignmentRSDB( std::ostream & output, 
		const Alignment * src );

/** enum describing the ways that two alignments can be combined
 * R: row
 * C: col
 */
typedef enum { RR, RC, CR, CC } CombinationMode;

/** return a pointer to an Alignment-object, that has been created by combining two others */
Alignment * combineAlignment( Alignment * dest, 
		const Alignment *src1, 
		const Alignment * src2, 
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
Alignment * copyAlignment( Alignment * dest, 
		const Alignment * src, 
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
Alignment * copyAlignmentRemoveRegion( Alignment * dest, 
		const Alignment * src, 
		Position row_from = NO_POS,
		Position row_to = NO_POS,
		Position col_from = NO_POS,
		Position col_to = NO_POS,
		Diagonal diagonal_from = 1,
		Diagonal diagonal_to = 0
);

/** copy one alignment into another using filter 
 */
Alignment * copyAlignment( Alignment * dest, 
		const Alignment * src, 
		const Alignment * filter, 
		const CombinationMode mode);

/** add one alignment to another
 */
Alignment * addAlignment2Alignment( Alignment * dest, const Alignment * src );

/** add one alignment to another. Map src using map_src2new.
 */
Alignment * addMappedAlignment2Alignment( Alignment * dest, 
		const Alignment * src, 
		const Alignment * map_src2new,
		const CombinationMode mode );

/** add one alignment to another. Map both row and column.
 */
Alignment * addMappedAlignments2Alignment( Alignment * dest, 
		const Alignment * src, 
		const Alignment * map_src_row2dest_row, 
		const Alignment * map_src_col2dest_col );


/** print a nice pairwise alignment */
void writePairAlignment( std::ostream & output, 
		const Alignandum * row, 
		const Alignandum * col, 
		const Alignment * ali );

/** write a nice pairwise alignment, allowing for wrapping around of the alignment */
void writeWraparoundAlignment( std::ostream & output, 
		const Alignandum * row, 
		const Alignandum * col, 
		const Alignment * ali,
		size_t max_insert_length = 30);

/** create an identity alignment between residues from and to in row using an offset for col */
Alignment * fillAlignmentIdentity( Alignment * dest, 
		Position row_from, 
		Position row_to, 
		Position col_offset = 0);

/** fill gaps in an alignment by doing a local alignment in each
     region.
 */
Alignment * fillAlignmentGaps( Alignment * dest,
		Alignator * alignator,
		Alignandum * row,
		Alignandum * col );

/** remove residues from an alignment, that are part of another alignment
     @param mode: specifies, which residues are looked up. If mode = RR, then every pairs is eliminated from dest,
     where the row is also present as a row-residue in filter.
 */
Alignment * filterAlignmentRemovePairs( Alignment * dest, 
		const Alignment * filter, 
		const CombinationMode mode );

/** remove residues from an alignment, that are part of another alignment
     @param mode: specifies, which residues are looked up. If mode = RR, then every pairs is eliminated from dest,
     where the row is also present as a row-residue in filter.
 */
Alignment * filterAlignmentRemovePairwiseSorted( Alignment * dest, 
		const Alignment * filter, 
		const CombinationMode mode );

/** fill a alignment given an explicit alignment 
     for example	row_ali = "AAA---CCCKKKAAA"
     col_ali = "AAAKKKCCCKKK--A"
     gaps are skipped, residue numbering starts at
     the first residue. Note that you will run into problems,
     if residues have been skipped in the explicit alignment.
 */
Alignment * fillAlignmentExplicit( Alignment * dest,
		Position row_from, 
		const std::string & row_ali,
		Position col_from, 
		const std::string & col_ali
);


/** fill a multiple alignment given compressed strings */
Alignment * fillAlignmentCompressed( Alignment * dest, 
		Position row_from, 
		const std::string & row_ali,
		Position col_from, 
		const std::string & col_ali
);


/** fill a multiple alignment given a compressed string in diagonal format. If 
     reverse is true, then row and column will be reversed while building alignment.
 */
Alignment * fillAlignmentCompressedDiagonal( Alignment * dest, 
		const std::string & ali,
		bool reverse = false
);


/** fill a multiple alignment from a stream. The format is a simple pairlist.
 */
Alignment * readAlignmentPairs( Alignment * dest, 
		std::istream & input,
		bool reverse = false
);


/** fill a multiple alignment given a compressed string in diagonal format but apply filter */
Alignment * fillAlignmentCompressedDiagonal( Alignment * dest, 
		const std::string & ali,
		Position row_from,
		Position row_to = NO_POS,
		Position col_from = NO_POS,
		Position col_to = NO_POS,
		Diagonal diagonal_from = -MAX_DIAGONAL,
		Diagonal diagonal_to = MAX_DIAGONAL
);

/** rescore alignment. This routine is generic as it uses the residue
     representation of row and col to calculate a score using a SubstitutionMatrix.
 */
Alignment * rescoreAlignment( Alignment * dest,
		const Alignandum * row,
		const Alignandum * col,
		const SubstitutionMatrix * matrix = NULL);


/** rescore alignment setting each pair to the same score 
 */
Alignment * rescoreAlignment( Alignment * dest,
		const Score score = 0);

/** rescore alignment. This routine tries to type-cast row and col and 
     accesses private data of row and col.
 */
Alignment * rescoreAlignmentPrivate( Alignment * dest,
		const Alignandum * row,
		const Alignandum * col,
		const SubstitutionMatrix * matrix = NULL);



/** calculate Alignment score given gap-penalties for row and column */
Alignment * calculateAffineScore( Alignment * dest, 
		Score gop, 
		Score gep );


/** fill an alignment with a repeat unit from a wrap-around alignment */
Alignment * fillAlignmentRepeatUnit( Alignment * dest, 
		const Alignment * source,
		Position first_row_residue = NO_POS,
		bool skip_negative_ends = false);


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

void fillAlignmentSummation( Alignment * dest1, 
		Alignment * dest2, 
		const Alignment * src,
		bool insert_gaps_row = true,
		bool insert_gaps_col = true, 
		bool use_end_row = false,
		bool use_end_col = false, 
		Position row_length = NO_POS,
		Position col_length = NO_POS);


/** complement a pairwise alignment. If there is a gap of the same length in both row and
     col, the corresponding residues are added to the alignment.
 */
Alignment * complementAlignment( Alignment * dest, 
		const Position max_length );

/** remove all those residues from an alignmnent, which are not
     in sequence. This ensures, that col_i < col_i+1 and row < row_i+1
     Only use with AlignmentVector

 */
Alignment * flattenAlignment( Alignment * dest );

/** split an alignment, if there are gaps larger than a certain threshold either in row or
     col or both.
 */
FragmentVector * splitAlignment( const Alignment * src, 
		const int max_gap_width,
		bool split_row = true,
		bool split_col = true);

/** split an alignment at points of intersection with another alignment.
 */ 
FragmentVector * splitAlignment( const Alignment * src1, 
		const Alignment * src2, 
		const CombinationMode mode );

/** starting from the ends of an alignment, remove 
    residues which do not contribute to a positive score.
 */
void pruneAlignment( Alignment * src,
		const Score gop,
		const Score gep);

/** calculate percent similarity of alignment */
double calculatePercentSimilarity( const Alignment * src);  

/** calculate percent identity of alignment. Since this depends on the objects mapped on the alignment, 
    you have to supply them. */
double calculatePercentIdentity (const Alignment * src, 
		const Alignandum * row, 
		const Alignandum * col); 


/** remove small fragments from alignment.
    This method removes fragments from an alignment. A fragment
    is a part of an alignment, that is short (max_fragment_length)
    and surrounded by large gaps (min_gap_length).
 */
void removeFragments( Alignment * dest,
		unsigned int window_length,
		unsigned int min_gap_length,
		Position row_length = NO_POS);

}

#endif	/* HELPERS_ALIGNATA_H */






