/*
  alignlib - a library for aligning protein sequences

  $Id: MultipleAlignment.h,v 1.6 2004/09/16 16:02:38 aheger Exp $

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

#ifndef MULTIPLEALIGNMENT_H
#define MULTIPLEALIGNMENT_H 1

#include <iosfwd>
#include <string>

#include "alignlib.h" 

namespace alignlib {


class Alignandum;
class Alignatum;
class Alignata;	
class Renderer;

/**
    @short Interface definition of mutltiple alignment objects.
    
    Multiple alignments are collection of aligned sequences (more specifically: objects of type
    @ref Alignata). A multiple alignment
    receives submitted aligned objects, or unaligned objects together with a multiple
    alignment, etc. Ownership belongs to the multiple alignment.
    
    Iterators for multiple alignments: not implemented yet, will come. So far, use
    operator[].
    
    The functions for obtaining an iterator over the multiple alignments is modelled
    after containers in the STL. I did not want to restrict the implementation to 
    vectors, (or lists, etc.), because some implementations might like to have it
    in a hash. But nevertheless, it should be possible to iterate over a multiple
    alignment and therefore we need a container.
    
    (in re maps: actually, it might be nice to then to templatize operator[](T&),
    so that access via id can be possible then as well. Before doing so, check for
    code bloat. On the other hand, when interacting via iterators, the whole class
    could be a template, i.e. MultipleAlignment< vector >, or MultipleAlignment< list >
    )
    
    The multiple alignment is a sort of container. Output-formatting, etc. is done
    by the @ref Alignata objects.
    
    The interface for this class is quite fat, because multiple alignments are used
    in a variety of contexts.

    This class is a protocol class and as such defines only the general interface

    About the parameter skip_gaps: When skip_gaps is set to true, then the current
    multiple alignment is regarded as a master multiple alignment and no gaps are 
    inserted into this alignment. When it is not set, gaps are added in the middle
    and the ends.

    10.4.2001: added default value NULL for parameter renderer in method registerRenderer() 

    @author Andreas Heger
    @version $Id: MultipleAlignment.h,v 1.6 2004/09/16 16:02:38 aheger Exp $
*/
class MultipleAlignment {

    friend std::ostream & operator<<( std::ostream &, const MultipleAlignment &);
    friend std::istream & operator>>( std::istream &, MultipleAlignment &);

    // class member functions
 public:
    /** iterator associated with multiple alignments */
    class iterator;
    
    /** const_iterator associated with multiple alignments */
    class const_iterator;

    // constructors and desctructors
    /** empty constructor */
    MultipleAlignment  ();

    /** copy constructor */
    MultipleAlignment  (const MultipleAlignment &);

    /** destructor */
    virtual ~MultipleAlignment ();

    //---------------------------------------------------------------------------------------
    /*------- accessors --------------------------------------------------------------------*/
    /** returns the length of the multiple alignment. All objects in a multiple alignment have
	the same length */
    virtual Position getLength() const = 0;

    /** sets the length of the multiple alignment. Raises an exception, if the mali is not empty
	*/
    virtual void setLength( Position length) = 0;

    /** returns the width of the multiple alignment, i.e., the number of objects in this multiple
	alignment */
    virtual int getWidth() const = 0;
    
    /** returns a const reference to the object at row in the multiple alignment. This allows treating
	the multiple alignment as a two-dimensional matrix. Since string does define operator[] as well, you
	can access the symbol in column c and row r by calling 
	symbol =  multiple_alignment[r][c]
    */
    virtual const std::string & operator[]( int row ) const = 0;		

    /** returns a pointer to the object at row from the multiple alignment. The pointer is not const,
	so you are free to do all sort of ugly stuff.*/
    virtual Alignatum * getRow( int row ) const = 0;

    /** erases an entry form the multiple alignment */
    virtual void eraseRow( int row ) = 0;

    /* ------------------ mutators ----------------------------------------------------------- */

    /*------------------- functions for adding new members to the multiple alignment---------*/
    

    /** add an aligned object to the multiple alignment. Ownership of this object is transferred to
	the multiple alignment. 
	@param src	       pointer to the aligned object to be added.
	@param alignment pointer to the alignment used for combining these two objects. If it is
		       not supplied, then it is assumed, that it is the identity alignment. In
		       that case src has to have the same length the multiple alignment. Note, the
		       multiple alignment is in col, the src is in row of the multiple alignment, so
		       when calling the member-function Alignata::Map() with a residue from
		       src, you get the correct position in the multiple alignment.
        @param mali_is_in_row		true, if the multiple alignment is in the row in alignment.
	@param insert_gaps_mali		true, if gaps shall be inserted into the multiple alignment.
	@param insert_gaps_alignatum	analogous to insert_gaps_alignatum.
	@param use_end_mali		true, if not-aligned residues at the ends of the multiple alignment
					shall be kept.
	@param use_end_alignatum	analogous to use_end_mali.
			      
    */
    virtual void add( Alignatum * src,
		      const Alignata * alignment = NULL,
		      bool mali_is_in_row = true,
		      bool insert_gaps_mali = true,
		      bool insert_gaps_alignatum= true,
		      bool use_end_mali = false,
		      bool use_end_alignatum = false) = 0;

    /** add the contents of a multiple alignment to the multiple alignment by mapping it through an alignment
     */
    virtual void add( const MultipleAlignment * src,
		      const Alignata * alignment = NULL,
		      bool mali_is_in_row = true,
		      bool insert_gaps_mali = true,
		      bool insert_gaps_alignatum= true,
		      bool use_end_mali = false,
		      bool use_end_alignatum = false) = 0;
    
    /** returns true, if there are no aligned objects in this alignment */
    virtual bool isEmpty() const = 0;

    /** clears the multiple alignment */
    virtual void clear() = 0;

    /** register a new renderer */
    virtual void registerRenderer( const Renderer * renderer = NULL) = 0;
    
    /** returns the consensus string for the multiple alignment */
    virtual std::string getConsensusString () const = 0;			

    /** returns a clone of this object */
    virtual MultipleAlignment * getClone() const = 0;
    
    /** returns an empty version of this object */
    virtual MultipleAlignment * getNew() const = 0;
    
    /** Write the multiple alignment to a stream

        @param segment_from
	@param segment_to
    
	Restrict output to a particular segment, if you wish so. Ranges
	start from 0 and are like Python-ranges, i.e., segment_to
	is not part of the range any more.
    */
    virtual void write( std::ostream & output, 
			Position segment_from = 0, 
			Position segment_to = 0) const = 0;

    /** Read a multiple alignment form a stream. I am not quite sure
	how this will work, because there might be different @ref Alignata
	objects in there.
     */
    virtual void read( std::istream & ) =0;
};

}

#endif /* MULTIPLE_ALIGNMENT_H */

