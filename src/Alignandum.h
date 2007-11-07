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
    Base class for objects that are to be aligned. Since all those objects are a sequence of 
    some sort, they have a length and fragments can be specified by two integer values. On request,
    they export const pointers to their member data for aligning.
    
    This class is responsible for mapping exporting only a segment of the sequence, if such is desired. 
    The alignment-algorithms always assume that an Alignandum object reaches from 1 to length. 
    
    This class is a protocol class and as such defines only the general interface

    @author Andreas Heger
    @version $Id: Alignandum.h,v 1.2 2004/01/07 14:35:31 aheger Exp $
    @short protocol class of alignable objects
*/


class Alignandum {
    /* friends ---------------------------------------------------------------------------- */
    friend  std::ostream & operator<<( std::ostream &, const Alignandum &);
    friend  std::istream & operator>>( std::istream &, Alignandum &);             
    
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

    /** get length of window */
    virtual Position	getLength() const = 0;

    /** use a segment for exporting and set segment to from and to 
        @param from     where segment starts
        @param to       where segment ends
    */
    virtual void useSegment( Position from = NO_POS, Position to = NO_POS) = 0;
    
    /** return first residue number in segment */
    virtual Position getFrom() const = 0;

    /** return last residue number in segment */
    virtual Position getTo() const = 0;

    /** return true if object is prepared for alignment (for cacheable types ) */
    virtual bool isPrepared() const = 0;	

    /** returns something, that can be used to access the member data. Fortunately it is now
	part of the standard, to change the return type of an inherited member function.
    */
    virtual const AlignandumData & getData() const = 0;
    
    /** get internal representation of residue in position pos */
    virtual Residue asResidue( Position pos ) const = 0;

    /** get internal representation of residue in position pos */
    virtual char asChar( Position pos ) const = 0;

    /** returns a string representation of the object */
    virtual std::string asString() const = 0;
    
    /** mask column at position x or in a segment y,
       if y is 0 
       */
    virtual void mask( Position from, Position to = 0) = 0;

    /** shuffle object */
    virtual void shuffle( unsigned int num_iteratinos = 1,
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

