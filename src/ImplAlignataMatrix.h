/*
  alignlib - a library for aligning protein sequences

  $Id: ImplAlignataMatrix.h,v 1.3 2004/03/19 18:23:40 aheger Exp $

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

#ifndef IMPL_ALIGNATA_MATRIX_H
#define IMPL_ALIGNATA_MATRIX_H 1

#include <iosfwd>
#include <vector>
#include "alignlib.h"
#include "ImplAlignata.h"

namespace alignlib 
{


class Alignandum;


//----------------------------------------------------------------------------------------------------------------------
/** @brief calculate the diagonal for a residue pair.
 */
inline Diagonal calculateDiagonal( const ResiduePAIR & p) { 
  return (p.mCol - p.mRow); 
}
//----------------------------------------------------------------------------------------------------------------------
/** @brief calculates the diagonal, so that (1,1) is (row_from, col_from)
 */
inline Diagonal calculateNormalizedDiagonal( const ResiduePAIR & p, 
						  Position row_from, 
						  Position col_from) { 
  return ( (p.mCol - col_from +1 ) - (p.mRow - row_from + 1)); 
}

/** @short base class for dotplots.

    Residues are kept in a list in memory. Note, that in this form of alignments,
    several column residues can be aligned to the same row residue and vice versa.
	
    Think of this Alignment as a dot-matrix. Actually, this is what it is used for
    in AlignatorDot. However, the full matrix is not stored memory, but as 
    a list of aligned residues (sorted first by row, then by column). The class provides
    and index for quick mapping of row to column (always the smallest col is returned).
    
    Note, that some of the functions meaningful in simple pairwise alignment object have
    a consistent but slightly different sense here. For example the alignment length
    returns the number of dots might therefore be longer than the sequences. 
    
    Also, mapRowToCol and mapColToRow give unpredictable results.
    
    Some implementation details:
    
    If the alignment is empty, mRowFrom is larger than mRowTo. This property is needed
    at some checks (for example in ImplAlignataMatrix_iterator)
    
    @author Andreas Heger
    @version $Id: ImplAlignataMatrix.h,v 1.3 2004/03/19 18:23:40 aheger Exp $
*/

class ImplAlignataMatrix : public ImplAlignata 
{

    friend class ConstIterator;
    friend class ImplAlignatorDots;

 public:

    typedef std::vector<ResiduePAIR *> PAIRVECTOR;
    typedef PAIRVECTOR::iterator PairIterator;
    typedef PAIRVECTOR::const_iterator PairConstIterator;

    //------------------> constructors / destructors <---------------------------------------------------------
    /** constructor */
    ImplAlignataMatrix( long ndots = 0);
    
    /** copy constructor */
    ImplAlignataMatrix( const ImplAlignataMatrix &src );

    /** destructor */
    virtual ~ImplAlignataMatrix();

    //------------------------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------------------------
    // this iterator iterates over all dots. 
    //!! to be implemented: another iterator, that iterates over rows
    // valid values for index are 0 to max_index - 1
    class ImplAlignataMatrix_ConstIterator : public Alignata::ConstIterator{
    public:
	
      ImplAlignataMatrix_ConstIterator( const PAIRVECTOR & pairs, 
					 Position index, 
					 Position max_index ) :
	mPairs( pairs ), mCurrentIndex(index), mMaximumIndex( max_index) {
      };
	
      ImplAlignataMatrix_ConstIterator( const ImplAlignataMatrix_ConstIterator & src ) : 
	mPairs( src.mPairs ), mCurrentIndex( src.mCurrentIndex ), mMaximumIndex( src.mMaximumIndex ) {};
      
      virtual ~ImplAlignataMatrix_ConstIterator() {};
      
      virtual ConstIterator * getClone() const {return new ImplAlignataMatrix_ConstIterator( *this );}
      
      /** dereference operator: runtime error, if out of bounds? */
      virtual const ResiduePAIR & getReference() const { 
	return (*mPairs[mCurrentIndex]);
      }
      
      /** for indirection */
      virtual const ResiduePAIR * getPointer() const { 
	  if (mCurrentIndex >= 0 && mCurrentIndex < mMaximumIndex)
	    return mPairs[mCurrentIndex]; 
	  else 
	    return NULL;
      }
      
      /** advance one position, until you find an aligned pair */
      virtual void next() { mCurrentIndex++; if (mCurrentIndex > mMaximumIndex) mCurrentIndex = mMaximumIndex; }
      
      /** step back one position, until you find an aligned pair */
      virtual void previous() { mCurrentIndex--; if (mCurrentIndex < -1) mCurrentIndex = -1; }
	
    private:
      const PAIRVECTOR & mPairs;
      Position mCurrentIndex;
      Position mMaximumIndex;
      
    };

    class ImplAlignataMatrix_Iterator : public Alignata::Iterator {
    public:
	
      ImplAlignataMatrix_Iterator( PAIRVECTOR & pairs, 
				   Position index, 
				   Position max_index ) :
	mPairs( pairs ), mCurrentIndex(index), mMaximumIndex( max_index) {
	};
	
	ImplAlignataMatrix_Iterator( const ImplAlignataMatrix_Iterator & src ) : 
	  mPairs( src.mPairs ), mCurrentIndex( src.mCurrentIndex ), mMaximumIndex( src.mMaximumIndex ) {};
	
	 virtual ~ImplAlignataMatrix_Iterator() {};
      
	virtual Iterator * getClone() const {return new ImplAlignataMatrix_Iterator( *this );}
	
	/** dereference operator: runtime error, if out of bounds? */
	virtual ResiduePAIR & getReference() const { 
	  return (*mPairs[mCurrentIndex]);
	}
 
	/** for indirection */
	virtual ResiduePAIR * getPointer() const { 
	  if (mCurrentIndex >= 0 && mCurrentIndex < mMaximumIndex)
	    return mPairs[mCurrentIndex]; 
	  else 
	    return NULL;
	}
      
	/** advance one position, until you find an aligned pair */
	virtual void next() { mCurrentIndex++; if (mCurrentIndex > mMaximumIndex) mCurrentIndex = mMaximumIndex; }
	
	/** step back one position, until you find an aligned pair */
	virtual void previous() { mCurrentIndex--; if (mCurrentIndex < -1) mCurrentIndex = -1; }
	
    private:
	PAIRVECTOR & mPairs;
	Position mCurrentIndex;
	Position mMaximumIndex;
	
    };
    
    /** return const iterator */
    virtual AlignataConstIterator begin() const; 
    
    /** return const iterator */
    virtual AlignataConstIterator end() const; 

    /** return const iterator */
    virtual AlignataIterator begin(); 
    
    /** return const iterator */
    virtual AlignataIterator end(); 

    //----------------> accessors <------------------------------------------------------------------------------

    /** returns the first residue aligned in row */
    virtual Position getRowFrom() const;

    /** returns the last residue aligned in row */
    virtual Position	getRowTo() const;

    /** returns the first residue aligned in col */
    virtual Position	getColFrom() const;

    /** returns the last residue aligned in col */
    virtual Position	getColTo()   const;

    /** returns the first aligned pair. Have to create a copy not a reference, because not all alignments will have
	a list of pairs ready */
    virtual ResiduePAIR front() const;
    
    /** returns the last aligned pair */
    virtual ResiduePAIR back() const;

    /** adds a pair of residue to the alignment */
    virtual void addPair( ResiduePAIR * new_pair ); 

    /** adds a pair of residues to the alignment */
    virtual void addPair( Position row, Position col, Score score = 0); 

    /** removes a pair of residues from the alignment */
    virtual void removePair( const ResiduePAIR & old_pair );

    /** retrieves a pair of residues from the alignment */
    virtual ResiduePAIR getPair( const ResiduePAIR & p) const;

    /** clear the current alignemnt */
    virtual void clear();

    /** move alignment */
    virtual void moveAlignment( Position row_offset, Position col_offset);

 protected:

    /** build index */
    virtual void buildIndex() const = 0;

    /** sort Dots by row and col */
    virtual void sortDots() const = 0;

    /** calculate alignment length */
    virtual void calculateLength() const;

    /** allocate Memory for ndots dots */
    virtual void allocateIndex( unsigned long size ) const;
    
    /** eliminate duplicates entries in mPairs. This only works, if mPairs has
	been sorted previously.
    */
    virtual void eliminateDuplicates() const;
    
    /** sort dots in range from-to by row */
    void sortDotsByRow(Position from, Position to) const;

    /** sort dots in range from-to by col */
    void sortDotsByCol(Position from, Position to) const;

    /** sort dots in range from from to to by diagonal */
    void sortDotsByDiagonal(Position from, Position to) const;

    /** List of residue pairs, mutatble, because they get sorted in-situ */
    mutable PAIRVECTOR mPairs;

    /** index of pairs for each row */
    mutable Dot * mIndex;

    /** coordinates of alignment */
    mutable Position mRowFrom;
    mutable Position mRowTo;
    mutable Position mColFrom;
    mutable Position mColTo;


 private:
    /* allocated size of size */
    mutable long mAllocatedIndexSize;

};

						  

}

#endif /* _ALIGNATA_H */

