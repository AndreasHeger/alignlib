/*
  alignlib - a library for aligning protein sequences

  $Id: ImplAlignatorFullDP.h,v 1.3 2004/03/19 18:23:41 aheger Exp $

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

#ifndef IMPL_ALIGNATOR_DP_H
#define IMPL_ALIGNATOR_DP_H 1

#include "alignlib.h"
#include "ImplAlignator.h"

namespace alignlib {

    class SubstitutionMatrix;
    class Alignandum;
    class Alignata;

/** codes in trace-back matrix */
#define TB_NOMATCH -999
#define TB_STOP       0
#define TB_INSERTION -1
#define TB_DELETION  -2
#define TB_MATCH     -3
#define TB_WRAP      -4


    /* re: Global functions and pointers for the fast determination of match score.
       
      I don't know how to use member functions as function pointers. After all, this is what
      inheritence is for. A possibility would be to automatically subclass an alignator-object.
      However, I do not like this idea, since this assumes that the parent has information 
      about the child. On the other hand, via the inlining mechanism a function call could be saved.
      The problem is when you want to change the algorithm by overloading align. Then the parent
      functions () do not know, what the child is. Therefore I use the static 
      functions. Maybe it is possible to separate the algorithm and the type-decision into different
      classes that interact.

      It should be possible, it is just a syntax problem?

      The danger of static functions is that the global pointers are unsafe, i.e. there exist just
      one copy for all alignator-objects, and my guess is that this code will never ever be threadsafe.
    */

/**
    @brief local,  dynamic programming alignment 
   
    This objects aligns two @ref Alignandum objects using a  dynamic programming algorithm with
    affine gap penalties. The objects to be aligned are always assumed to be starting at residue 1.
    The object @ref Alignandum are responsible for mapping windows back and forth.
    
    This class implements the back-tracking algorithm used by several of its children.

    @author Andreas Heger
    @version $Id: ImplAlignatorDP.h,v 1.3 2004/03/19 18:23:41 aheger Exp $
*/
class ImplAlignatorDP : public ImplAlignator {

    typedef Score (ImplAlignatorDP::*MATCH_FUNCTION_POINTER)( Position, Position); 
		      
 public:
    /* Constructors and destructors */

    /** set affine gap penalties and substitution matrix 
     @param subst_matrix	pointer to substitution matrix
     @param row_gop		gap opening penalty in row
     @param row_gep		gap elongation penalty in row
     @param col_gop		gap opening penalty in column, default = row
     @param col_gep		gap elongation penalty in row, default = col
     
    */
    ImplAlignatorDP(  const SubstitutionMatrix * subst_matrix ,
			  Score row_gop, Score row_gep, 
			  Score col_gop = 0,Score col_gep = 0 );
			

    /** copy constructor */
    ImplAlignatorDP( const ImplAlignatorDP & );

    /** destructor */
    virtual ~ImplAlignatorDP();

    /* operators------------------------------------------------------------------------------ */
    /** method for aligning two arbitrary objects */
    virtual Alignata * align(const Alignandum *, const Alignandum *, Alignata *);

    /** align two objects and just return the score */
    virtual Score getAlignmentScore(const Alignandum *, const Alignandum * );		

    /* member access functions--------------------------------------------------------------- */

 protected:
    /** perform initialization before alignment */
    virtual void startUp(const Alignandum * row, const Alignandum * col, Alignata * ali);                     
    
    /** clean up temporary memory after alignment step */
    virtual void cleanUp(const Alignandum * row, const Alignandum * col, Alignata * ali);                     

    /** traces back through trace matrix and put in the alignment in Alignata-object */
    virtual void traceBack( const Alignandum * row, const Alignandum * col, Alignata * result);				

    /** return index for given row and length */
    inline int getTraceIndex( int r, int c ) { return (c * mRowLength + r); }; 

    /** perform the alignment */
    virtual void performAlignment(const Alignandum * row, const Alignandum *col, Alignata * result);

    /** perform the alignment */
    virtual void performAlignmentWithoutTraceBack(const Alignandum * row, const Alignandum *col);
 private:
    /** Match-functions */
    Score matchSequenceSequence( Position row, Position col );
    Score matchSequenceSequenceTCO( Position row, Position col );
    Score matchSequenceTCOSequence( Position row, Position col );
    Score matchProfileSequence( Position row, Position col );
    Score matchSequenceProfile( Position row, Position col );
    Score matchProfileProfile( Position row, Position col );

    /* member data --------------------------------------------------------------------------- */
 protected:
    /** pointer to the trace matrix */
    int *mTrace;

    /** row, where trace ended */
    int mRowLast;

    /** column, where trace ended */
    int mColLast;		       

    /** maximum score in matrix */
    Score mScore;

    /** internal helper array for the calculation of affine gap penalties */
    Score *mCC;			

    /** internal helper array for the calculation of affine gap penalties */
    Score *mDD;
    
    /** pointer to function that calculates match-score */
    MATCH_FUNCTION_POINTER mMatchFunction;

    // pointers to memory location of encoded sequences/profiles/...

    /** pointer to member data of row/col : AlignandumSequence */
    const Residue * mRowSequence; 
    const Residue * mColSequence;
    
    /** pointer to member data of row/col AlignandumProfile */
    const ProfileColumn * mRowProfile;
    const ProfileColumn * mColProfile;
    const FrequencyColumn * mRowFrequencies;
    const FrequencyColumn * mColFrequencies;
    
    /** pointer to tco-columns, only row, because subst-matrix is not symmetric */
    const TCOMATRIXCOLUMN * mRowSequenceTCO;

    /** length of row plus 1, saved here for quick calculation in getTraceIndex */
    Position mRowLengthPlus1;
};


}

#endif /* _ALIGNATOR_MATRIX_H */

