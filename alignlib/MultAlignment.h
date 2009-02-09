/*
  alignlib - a library for aligning protein sequences

  $Id: MultAlignment.h,v 1.6 2004/09/16 16:02:38 aheger Exp $

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

#ifndef MULTALIGNMENT_H
#define MULTALIGNMENT_H 1

#include <iosfwd>
#include <string>

#include "alignlib_fwd.h"

namespace alignlib
{

/**
    @short Protocoll class for multiple alignments.

    Multiple alignments are collections of aligned sequences (objects of type
    @ref Alignatum).

    This class is a protocol class and as such defines only the general interface

    About the parameter skip_gaps: When skip_gaps is set to true, then the current
    multiple alignment is regarded as a master multiple alignment and no gaps are
    inserted into this alignment. When it is not set, gaps are added in the middle
    and the ends.

    @author Andreas Heger
    @version $Id: MultAlignment.h,v 1.6 2004/09/16 16:02:38 aheger Exp $
*/
class MultAlignment
{

    friend std::ostream & operator<<( std::ostream &, const MultAlignment &);

    // class member functions
 public:

    // constructors and desctructors
    /** empty constructor */
    MultAlignment  ();

    /** copy constructor */
    MultAlignment  (const MultAlignment &);

    /** destructor */
    virtual ~MultAlignment ();

    //---------------------------------------------------------------------------------------
    /*------- accessors --------------------------------------------------------------------*/
    /** returns the length (number of columns) of the multiple alignment.
     *
     * All objects in a multiple alignment have the same length.
     *
     * @return the length (aligned positions) of the multiple alignment. */
    virtual Position getLength() const = 0;

    /** returns the number of sequences in this multiple alignment.
     *
     * @return number of sequences in alignment.
     */
    virtual int getNumSequences() const = 0;

    /** returns a row of the multiple alignment.
     *
     * @param row row of multiple alignment.
     * return multiple alignment
    */
    virtual const HAlignment operator[]( int row ) const = 0;

    /** returns a row of the multiple alignment.
     *
     * @param row row of multiple alignment.
     * return a @ref Alignatum object
    */
    virtual const HAlignment getRow( int row ) const = 0;

    /** erases an row from the multiple alignment
     *
     * @param row row to erase.
     * */
    virtual void eraseRow( int row ) = 0;

    /** return true, if a column is aligned.
     *
     * Unaligned columns result from adding
     * new sequences to the multiple alignment.
     *
     * @return true, if column @col is aligned.
     * */
    virtual bool isAligned( const Position & col ) = 0;


    /* ------------------ mutators ----------------------------------------------------------- */

    /** expand multiple alignment
     *
     * This will add columns into the multiple alignment for unaligned positions
     * in the rows.
     *
     * @param sequences  	a list of sequences. If this array is not null,
     * 						sequence lengths from this array will be used to expand
     * 						the multiple alignment before the first and after the last
     * 						column.
     * */
    virtual void expand(
    		const HAlignandumVector & sequences ) = 0;

    /*------------------- functions for adding new members to the multiple alignment---------*/


    /** add an @ref Alignment object to the multiple alignment.
     *
     * The alignment object maps the sequence to multiple alignment columns.
     *
	@param src	 @ref Alignment object to add.
	@param alignment @ref Alignment that maps src to mali.
    */
    virtual void add(
    		const HAlignment & map_mali2sequence ) = 0;

    /** add a @ref MultAlignment object to the multiple alignment.
     *
     * The alignment object maps the sequence to multiple alignment columns.
     * Note that some alignment information can be potentially lost. If two
     * sequence positions are aligned in @param src, but that column is not
     * in map_mali2sequnce, then the alignment of these two residues is lost.
     *
	@param src	 @ref Alignment object to add.
	@param alignment @ref Alignment that maps src to mali.
    */
    virtual void add(
    		const HMultAlignment & src,
    		const HAlignment & map_mali2sequence ) = 0;

    /** returns true, if the alignment is empty.
     *
     * @return true, if the alignment is emtpy.
     * */
    virtual bool isEmpty() const = 0;

    /** clears the multiple alignment
    */
    virtual void clear() = 0;

    /** returns a clone of this object
     * @return a copy of this object.
     * */
    virtual HMultAlignment getClone() const = 0;

    /** returns a new object of this type.
     * @return a new object of this type.
     * */
    virtual HMultAlignment getNew() const = 0;

    /** write the multiple alignment to a stream
     *
     * @param output output stream.
    */
    virtual void write( std::ostream & output ) const = 0;
};

}

#endif /* MULTIPLE_ALIGNMENT_H */

