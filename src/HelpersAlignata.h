/*
  alignlib - a library for aligning protein sequences

  $Id: HelpersAlignata.h,v 1.3 2004/10/14 23:34:09 aheger Exp $

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
#include <vector>
#include "alignlib.h"

namespace alignlib {
    
    /** Helper functions for class Alignata:
	
	1. factory functions
	
	2. accessor functions for default objects
	
	3. convenience functions
    */

    class Alignator;
    class Alignata;
    class Alignandum;
    class SubstitutionMatrix;

    typedef std::vector<Alignata*> FragmentVector;
    
 /* -------------------------------------------------------------------------------------------------------------------- */
 /* 1. factory functions */
    /** return a pointer to an Alignata-object */
    Alignata * makeAlignataSet(); 

    /** return a pointer to an Alignata-object (iterate column-wise)*/
    Alignata * makeAlignataSetCol(); 

    /** return a pointer to an Alignata-object */
    Alignata * makeAlignataHash(); 

    /** return a pointer to an Alignata-object */
    Alignata * makeAlignataHashDiagonal(); 
    
    /** return a pointer to an Alignata-object */
    Alignata * makeAlignataVector(); 

    /** return a pointer to a new Alignata-object */
    Alignata * makeAlignataMatrixRow( long ndots = 0);

    /** return a pointer to a new Alignata-object */
    Alignata * makeAlignataMatrixUnsorted( long ndots = 0);

    /** return a pointer to a new Alignata-object */
    Alignata * makeAlignataMatrixDiagonal( long ndots = 0);

    
 /* -------------------------------------------------------------------------------------------------------------------- */
 /* 2. accessor functions for default objects */

    
 /* -------------------------------------------------------------------------------------------------------------------- */
 /* 3. convenience functions */


    /** Print aligment in table format */
    void writeAlignataTable( std::ostream & output, 
			     const Alignata * src,
			     unsigned int ncols = 8,
			     bool with_scores = true);

 /** write compressed alignment into streams. The alignment can be restricted to
     a region specifying the columns. If you want to restrict it to row-residues,
     you can get those boundaries by calling mapRowToCol

     @param src	Alignment
     @param out_row	string in which to store the row in compressed format
     @param out_col	string in which to store the col in compressed format
     @param col_from	beginning of segment to use
     @param col_to	end of segment to use
  */
    void writeAlignataCompressed( const Alignata * src, 
				  std::string & out_row , 
				  std::string & out_col,
				  Position col_from = 0,
				  Position col_to = 0
				  );

 /** write compressed alignment into streams. The alignment can be restricted to
     a region specifying the columns. A further filter can be
     applied, that only saves a tube. 

     If diagonal_from is larger than diagonal_to, then the whole range is used.

     @param src	Alignment
     @param out_row	string in which to store the row in compressed format
     @param reverse  whether to reverse row and column
     @param row_from	beginning of segment to use
     @param row_to	end of segment to use
     @param col_from	beginning of segment to use
     @param col_to	end of segment to use
     @param diagonal_form beginning of tube to use
     @param diagonal_to end of tube to use
  */
 void writeAlignataCompressedDiagonal( const Alignata * src, 
				       std::string & out_row, 
				       bool reverse = false,
				       Position row_from = 0,
				       Position row_to = 0,
				       Position col_from = 0,
				       Position col_to = 0,
				       Diagonal diagonal_from = MAX_DIAGONAL,
				       Diagonal diagonal_to = -MAX_DIAGONAL
				       );

 /** write an alignment in rsdb-format */
 void writeAlignataRSDB( std::ostream & output, const Alignata * src );

 /** return a pointer to an Alignata-object, that has been created by combining two others */
 typedef enum { RR, RC, CR, CC } COMBINATION_MODE;

 Alignata * combineAlignata( Alignata * dest, const Alignata *src1, const Alignata * src2, const COMBINATION_MODE mode);

 /** map row/column of an alignment to the combination of another alignment */
 Alignata * combineAlignataMali( Alignata * map_ali2mali, 
				 const Alignata * map_ali2row, 
				 const Alignata * map_row2col, 
				 const COMBINATION_MODE mode);
 
 /** copy one alignment into another. I can not do this using a copy constructor, since virtual functions are not
     resolved by then. I also do not want a clone, for example if I want to change the underlying implementation

     @param row_from	beginning of segment to use
     @param row_to	end of segment to use
     @param col_from	beginning of segment to use
     @param col_to	end of segment to use
     @param diagonal_form beginning of tube to use
     @param diagonal_to end of tube to use
 */
 Alignata * copyAlignata( Alignata * dest, 
			  const Alignata * src, 
			  Position row_from = 0,
			  Position row_to = 0,
			  Position col_from = 0,
			  Position col_to = 0,
			  Diagonal diagonal_from = -MAX_DIAGONAL,
			  Diagonal diagonal_to = MAX_DIAGONAL
			  );

 /** copy one alignment into another. In contrast to copyAlignata, this function keeps the dots outside
     of region.
     @param row_from	beginning of segment to use
     @param row_to	end of segment to use
     @param col_from	beginning of segment to use
     @param col_to	end of segment to use
     @param diagonal_form beginning of tube to use
     @param diagonal_to end of tube to use
 */
 Alignata * copyAlignataRemoveRegion( Alignata * dest, 
				      const Alignata * src, 
				      Position row_from = 0,
				      Position row_to = 0,
				      Position col_from = 0,
				      Position col_to = 0,
				      Diagonal diagonal_from = 1,
				      Diagonal diagonal_to = 0
				      );

 /** copy one alignment into another using filter 
  */
  Alignata * copyAlignata( Alignata * dest, 
			   const Alignata * src, 
			   const Alignata * filter, 
			   const COMBINATION_MODE mode);

 /** add one alignment to another
  */
 Alignata * addAlignata2Alignata( Alignata * dest, const Alignata * src );

 /** add one alignment to another. Map src using map_src2new.
  */
 Alignata * addMappedAlignata2Alignata( Alignata * dest, 
					const Alignata * src, 
					const Alignata * map_src2new,
					const COMBINATION_MODE mode );

 /** add one alignment to another. Map both row and column.
  */
 Alignata * addMappedAlignatas2Alignata( Alignata * dest, 
					 const Alignata * src, 
					 const Alignata * map_src_row2dest_row, 
					 const Alignata * map_src_col2dest_col );
					 

 /** print a nice pairwise alignment */
 void writePairAlignment( std::ostream & output, 
			  const Alignandum * row, 
			  const Alignandum * col, 
			  const Alignata * ali );

 /** write a nice pairwise alignment, allowing for wrapping around of the alignment */
 void writeWraparoundAlignment( std::ostream & output, 
				const Alignandum * row, 
				const Alignandum * col, 
				const Alignata * ali );

 /** create an identity alignment between residues from and to in row using an offset for col */
 Alignata * fillAlignataIdentity( Alignata * dest, 
				  Position row_from, 
				  Position row_to, 
				  Position col_offset = 0);
 
 /** fill gaps in an alignment by doing a local alignment in each
     region.
 */
 Alignata * fillAlignataGaps( Alignata * dest,
			      Alignator * alignator,
			      Alignandum * row,
			      Alignandum * col );

 /** remove residues from an alignment, that are part of another alignment
     @param mode: specifies, which residues are looked up. If mode = RR, then every pairs is eliminated from dest,
     where the row is also present as a row-residue in filter.
 */
 Alignata * filterAlignataRemovePairsCopy( Alignata * dest, 
					   const Alignata * filter, 
					   const COMBINATION_MODE mode );

 /** remove residues from an alignment, that are part of another alignment
     @param mode: specifies, which residues are looked up. If mode = RR, then every pairs is eliminated from dest,
     where the row is also present as a row-residue in filter.
 */
 Alignata * filterAlignataRemovePairs( Alignata * dest, 
				       const Alignata * filter, 
				       const COMBINATION_MODE mode );

 /** remove residues from an alignment, that are part of another alignment
     @param mode: specifies, which residues are looked up. If mode = RR, then every pairs is eliminated from dest,
     where the row is also present as a row-residue in filter.
 */
 Alignata * filterAlignataRemovePairwiseSorted( Alignata * dest, 
						const Alignata * filter, 
						const COMBINATION_MODE mode );
 
 /** fill a alignment given an explicit alignment 
     for example	row_ali = "AAA---CCCKKKAAA"
     col_ali = "AAAKKKCCCKKK--A"
     gaps are skipped, residue numbering starts at
     the first residue. Note that you will run into problems,
     if residues have been skipped in the explicit alignment.
 */
 Alignata * fillAlignataExplicit( Alignata * dest,
				  Position row_from, 
				  const std::string & row_ali,
				  Position col_from, 
				  const std::string & col_ali
				  );
 

 /** fill a multiple alignment given compressed strings */
 Alignata * fillAlignataCompressed( Alignata * dest, 
				    Position row_from, 
				    const std::string & row_ali,
				    Position col_from, 
				    const std::string & col_ali
				    );


 /** fill a multiple alignment given a compressed string in diagonal format. If 
     reverse is true, then row and column will be reversed while building alignment.
  */
 Alignata * fillAlignataCompressedDiagonal( Alignata * dest, 
					    const std::string & ali,
					    bool reverse = false
					    );


 /** fill a multiple alignment from a stream. The format is a simple pairlist.
  */
 Alignata * readAlignataPairs( Alignata * dest, 
			       std::istream & input,
			       bool reverse = false
			       );


 /** fill a multiple alignment given a compressed string in diagonal format but apply filter */
 Alignata * fillAlignataCompressedDiagonal( Alignata * dest, 
					    const std::string & ali,
					    Position row_from,
					    Position row_to = 0,
					    Position col_from = 0,
					    Position col_to = 0,
					    Diagonal diagonal_from = -MAX_DIAGONAL,
					    Diagonal diagonal_to = MAX_DIAGONAL
					    );

 /** rescore alignment. This routine is generic as it uses the residue
     representation of row and col to calculate a score using a SubstitutionMatrix.
  */
 Alignata * rescoreAlignment( Alignata * dest,
			      const Alignandum * row,
			      const Alignandum * col,
			      const SubstitutionMatrix * matrix = NULL);


/** rescore alignment setting each pair to the same score 
  */
 Alignata * rescoreAlignment( Alignata * dest,
			      const Score score = 0);

 /** rescore alignment. This routine tries to type-cast row and col and 
     accesses private data of row and col.
  */
 Alignata * rescoreAlignmentPrivate( Alignata * dest,
				     const Alignandum * row,
				     const Alignandum * col,
				     const SubstitutionMatrix * matrix = NULL);



 /** calculate Alignment score given gap-penalties for row and column */
 Alignata * calculateAffineScore( Alignata * dest, 
				  Score gop, 
				  Score gep );
 

 /** fill an alignment with a repeat unit from a wrap-around alignment */
 Alignata * fillAlignataRepeatUnit( Alignata * dest, 
				    const Alignata * source,
				    Position first_row_residue = 0,
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

void fillAlignataSummation( Alignata * dest1, 
			    Alignata * dest2, 
			    const Alignata * src,
			    bool insert_gaps_row = true,
			    bool insert_gaps_col = true, 
			    bool use_end_row = false,
			    bool use_end_col = false, 
			    Position row_length = 0,
			    Position col_length = 0);


 /** complement a pairwise alignment. If there is a gap of the same length in both row and
     col, the corresponding residues are added to the alignment.
 */
 Alignata * complementAlignata( Alignata * dest, 
				const Position max_length );
 
 /** remove all those residues from an alignmnent, which are not
     in sequence. This ensures, that col_i < col_i+1 and row < row_i+1
     Only use with AlignataVector

 */
Alignata * flattenAlignata( Alignata * dest );

 /** split an alignment, if there are gaps larger than a certain threshold either in row or
     col or both.
 */
FragmentVector * splitAlignata( const Alignata * src, 
				const int max_gap_width,
				bool split_row = true,
				bool split_col = true);

/** split an alignment at points of intersection with another alignment.
 */ 
FragmentVector * splitAlignata( const Alignata * src1, 
				const Alignata * src2, 
				const COMBINATION_MODE mode );

/** starting from the ends of an alignment, remove 
    residues which do not contribute to a positive score.
*/
void pruneAlignata( Alignata * src,
		    const Score gop,
		    const Score gep);
 
/** calculate percent similarity of alignment */
double calculatePercentSimilarity( const Alignata * src);  
    
/** calculate percent identity of alignment. Since this depends on the objects mapped on the alignment, 
    you have to supply them. */
double calculatePercentIdentity (const Alignata * src, 
				 const Alignandum * row, 
				 const Alignandum * col); 


/** remove small fragments from alignment.
    This method removes fragments from an alignment. A fragment
    is a part of an alignment, that is short (max_fragment_length)
    and surrounded by large gaps (min_gap_length).
 */
void removeFragments( Alignata * dest,
		      unsigned int window_length,
		      unsigned int min_gap_length,
		      Position row_length = 0);

}

#endif	/* HELPERS_ALIGNATA_H */






