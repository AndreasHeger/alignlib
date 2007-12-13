/*
  alignlib - a library for aligning protein sequences

  $Id: ImplMultipleAlignmentDots.h,v 1.2 2004/01/27 12:14:49 aheger Exp $

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

#ifndef IMPL_MULTIPLE_ALIGNMENT_DOTS_H
#define IMPL_MULTIPLE_ALIGNMENT_DOTS_H 1

#include <iosfwd>
#include <string>
#include <list>

#include "alignlib.h" 
#include "MultipleAlignment.h"

namespace alignlib 
{

/** 
    Multiple alignments are collection of aligned sequences (more specifically: objects of type
    @ref Alignment). A multiple alignment
    receives submitted aligned objects, or unaligned objects together with a multiple
    alignment, etc. Ownership belongs to the multiple alignment.

    The multiple alignment is a sort of container. Output-formatting, etc. is done
    by the @ref Alignment-objects.
    
    The interface for this class is quite fat, because multiple alignments are used
    in a variety of contexts.

    In this implementation the rows are stored as pairs of sequences and alignments.
    Adding alignments is quick.
    
    @author Andreas Heger
    @version $Id: ImplMultipleAlignmentDots.h,v 1.2 2004/01/27 12:14:49 aheger Exp $
    @short an implementation for multiple alignments
*/

class Alignandum;
class Alignatum;
class Alignment;	
class Renderer;

/** This structure takes possession of the passed objects. 
 */
struct MaliRow 
{
  MaliRow(); 

  MaliRow( Alignatum * input, Alignment * map_alignatum2mali, Alignatum * output = NULL);

  ~MaliRow();
  
  /** the sequence */
  Alignatum * mAlignatumInput; 
  /** the dots */
  Alignment * mMapMali2Alignatum;
  /** the rendered Alignatum. Needed, for returning as reference */
  Alignatum * mAlignatumOutput;
};

class ImplMultipleAlignmentDots : public MultipleAlignment 
{

  // class member functions
 public:

    // constructors and desctructors
    /** empty constructor 
     */
    ImplMultipleAlignmentDots  (const bool compress_unaligned_columns = true,
				const int max_insertion_length = -1);

    /** copy constructor */
    ImplMultipleAlignmentDots  (const ImplMultipleAlignmentDots &);

    /** destructor */
    virtual ~ImplMultipleAlignmentDots ();

    //---------------------------------------------------------------------------------------
    /*------- accessors --------------------------------------------------------------------*/
    /** returns the length of the multiple alignment. All objects in a multiple alignment have
	the same length */
    virtual Position getLength() const;

    /** sets the length of the multiple alignment. Raises an exception, if the mali is not empty
	*/
    virtual void setLength( Position length);

    /** returns the width of the multiple alignment, i.e., the number of objects in this multiple
	alignment */
    virtual int getWidth() const;
    
    /** returns a const reference to the object at row in the multiple alignment. This allows treating
	the multiple alignment as a two-dimensional matrix. Since string does define operator[] as well, you
	can access the symbol in column c and row r by calling 
	symbol =  multiple_alignment[r][c]
    */
    virtual const std::string & operator[]( int row ) const;		

    /** returns a pointer to the object at row from the multiple alignment. The pointer is not const,
	so you are free to do all sort of ugly stuff.*/
    virtual Alignatum * getRow( int row ) const;
    
    /** erases an entry form the multiple alignment */
    virtual void eraseRow( int row );

    /* ------------------ mutators ----------------------------------------------------------- */

    /*------------------- functions for adding new members to the multiple alignment---------*/

    /** add an aligned object to the multiple alignment. Ownership of this object is transferred to
	the multiple alignment. 
	@param src	       pointer to the aligned object to be added.
	@param alignment pointer to the alignment used for combining these two objects. If it is
		       not supplied, then it is assumed, that it is the identity alignment. In
		       that case src has to have the same length the multiple alignment. Note, the
		       multiple alignment is in col, the src is in row of the multiple alignment, so
		       when calling the member-function Alignment::Map() with a residue from
		       src, you get the correct position in the multiple alignment.
	@param skip_gaps false, if gaps in existing multiple alignment shall be introduced
    */
    virtual void add( Alignatum * src,
		      const Alignment * alignment = NULL,
		      bool mali_is_in_row = true,
		      bool insert_gaps_mali = true,
		      bool insert_gaps_alignatum= true,
		      bool use_end_mali = false,
		      bool use_end_alignatum = false);

    /** add the contents of a multiple alignment to the multiple alignment by mapping it through an alignment
	@param src	       pointer to the unaligned object to be added.
	@param alignment pointer to the alignment used for combining these two objects. If it is
		       not supplied, then it is assumed, that it is the identity alignment. In
		       that case src has to have the same length the multiple alignment. Note, the
		       multiple alignment is in col, the src is in row of the multiple alignment, so
		       when calling the member-function Alignment::Map() with a residue from
		       src, you get the correct position in the multiple alignment.
          @param skip_gaps false, if gaps in existing multiple alignment shall be introduced
			   
     */
    virtual void add( const MultipleAlignment * src,
		      const Alignment * alignment = NULL,
		      bool mali_is_in_row = true,
		      bool insert_gaps_mali = true,
		      bool insert_gaps_alignatum= true,
		      bool use_end_mali = false,
		      bool use_end_alignatum = false);
    
    /** returns the consensus string for the multiple alignment */
    virtual std::string getConsensusString () const ;			

    /** register a new renderer */
    virtual void registerRenderer( const Renderer * renderer );

    /** returns true, if there are no aligned objects in this alignment */
    virtual bool isEmpty() const;

    /** clears the multiple alignment */
    virtual void clear();

    /** returns a clone of this object */
    virtual MultipleAlignment * getClone() const;
    
    /** returns an empty version of this object */
    virtual MultipleAlignment * getNew() const;

    /** Write the multiple alignment to a stream
     */
    virtual void write( std::ostream & output, 
			Position segment_from = 0, 
			Position segment_to = 0) const ;
    
 protected:
    /** render the multiple alignment */
    virtual void updateRows() const;
    
    /** free all memory. Tell all stored objects to destruct themselves */
    virtual void freeMemory();

 private:
    /** the length of the multiple alignment */
    mutable int mLength;                       
    
    /** I store an array of vectors. The pointers can not be const, because sequences are told to rescale. */
    mutable std::vector< MaliRow > mRows;             

    /** the Renderer */
    const Renderer * mRenderer;

    /** whether or not to compress unaligned columns */
    bool mCompressUnalignedColumns;

    /** maximal length of an insertion */
    int mMaxInsertionLength;
    
};


}

#endif /* IMPL_MULTIPLE_ALIGNMENT_DOTS_H */

