/*
  alignlib - a library for aligning protein sequences

  $Id: Alignandum.h,v 1.2 2004/01/07 14:35:31 aheger Exp $

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

#ifndef ALIGNANDUM_H
#define ALIGNANDUM_H 1

#include <iosfwd>
#include <string>
#include "alignlib.h"

namespace alignlib {

/** this is an empty class definition, that will be subclassed by other descendents
    of Alignandum. I have to do this, so that there will be no warning messages. This
    class contains handles (const pointers) to the member data of Alignandum objects.
*/

struct AlignandumData{};

/** 
    Base class for objects that are to be aligned, typically sequences or profiles.
        
    The class can restrict access to a sequence to a sequene range. Ranges are given
    in open/closed notation starting from 0. Thus, the segment 0..5 includes residues
    0,1,2,3,4.
    
    This class is a protocol class and as such defines only the general interface

    On request, they export const pointers to their member data for aligning. This
    is taken care of behind the scences.
    
    @author Andreas Heger
    @version $Id: Alignandum.h,v 1.2 2004/01/07 14:35:31 aheger Exp $
    @short protocol class of alignable objects
*/


class Alignandum 
{
    /* friends ---------------------------------------------------------------------------- */
    friend  std::ostream & operator<<( std::ostream &, const Alignandum &);
    
 public:
    /* constructors and desctructors------------------------------------------------------- */
    
    /** empty constructor */
    Alignandum();
    
    /** copy constructor */
    Alignandum( const Alignandum &);

    /** desctructor */
    virtual ~Alignandum();

    /** accessors ------------------------------------------------------------------------- */

    /** return an identical copy of this object */
    virtual Alignandum * getClone() const = 0;

    /** get length of sequence */
    virtual Position	getLength() const = 0;

    /** restrict the use of the sequence to a segment. If no coordinates are given,
     * the full sequence is used.
     *  
        @param from     where segment starts
        @param to       where segment ends
    */
    virtual void useSegment( Position from = NO_POS, Position to = NO_POS) = 0;
    
    /** return first residue number in segment */
    virtual Position getFrom() const = 0;

    /** return last residue number in segment */
    virtual Position getTo() const = 0;

    /** return true if object is prepared for alignment (for cacheable types ) 
     * This function permits lazy evaluation of of some alignable types like 
     * profiles. */
    virtual bool isPrepared() const = 0;	

    /** returns a structure that can be used to access internal data.
    */
    virtual const AlignandumData & getData() const = 0;
    
    /** get internal representation of residue in position pos 
     */
    virtual Residue asResidue( Position pos ) const = 0;

    /** get character representation of residue in position pos
     * using the default translator. 
     */
    virtual char asChar( Position pos ) const = 0;

    /** returns a string representation of the object 
     */
    virtual std::string asString() const = 0;
    
    /** mask column at position x or, if y is not omitted,
     * in a range.
       */
    virtual void mask( Position from, Position to = NO_POS) = 0;

    /** shuffle object 
     * */
    virtual void shuffle( unsigned int num_iterations = 1,
			  Position window_size = 0 ) = 0;

    /* Mutators ------------------------------------------------------------------------------ */

    /** load data into cache, if cacheable type */
    virtual void prepare() const = 0;						
    
    /** discard cache, if cacheable type */
    virtual void release() const = 0;					       

    /** write human readable output to stream.
     */
    virtual void write( std::ostream & output ) const = 0;

    /** save state of object into stream
     */
    virtual void save( std::ostream & input ) const = 0;
        
};

}
#endif /* ALIGNANDUM_H */

