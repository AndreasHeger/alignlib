/*
  alignlib - a library for aligning protein sequences

  $Id: Alignatum.h,v 1.5 2004/03/19 18:23:39 aheger Exp $

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

#ifndef ALIGNATUM_H
#define ALIGNATUM_H 1

#include <iosfwd>
#include <string>
#include "alignlib_fwd.h"
#include "alignlib_fwd.h"

namespace alignlib 
{

/** @short interface definition of Alignatum objects.

    Objects of this type represent an aligned object. This class is responsible for
    producing a nice, real-world representation of the aligned object. To this end, 
    this class provides several methods for modifying an aligned string.

    There is no method for inserting residues here. Do this, by modifying alignments
    before creating an Alignatum-object.
    
    @author Andreas Heger
    @version $Id: Alignatum.h,v 1.5 2004/03/19 18:23:39 aheger Exp $
*/
class Alignatum 
{
    // class member functions
    friend std::ostream & operator<<( std::ostream &, const Alignatum &);
    
 public:
    // constructors and desctructors
    /** when building a new Alignatum sequence by cloning and an alignment, this specifies,
	if the new sequence is in row or column
    */

    /** constructor */
    Alignatum  ();

    /** copy constructor */
    Alignatum  (const Alignatum &);

    /* destructor */
    virtual ~Alignatum ();

    /*-----> accessors <----------------------------------------------------- */

    /** write object into stream nicely formatted */
    virtual void writeRow( 
    		std::ostream & output, 
    		const HRenderer & renderer,
    		Position segment_start = 0, 
    		Position segment_end = 0 ) const = 0;

    /** return a copy of the string representation */
    virtual std::string getString() const = 0;

    /** return the string representation */
    virtual const std::string & getStringReference() const = 0;
    
    /** return the number of the first residue */
    virtual Position getFrom() const = 0;

    /** return the number of the last residue */
    virtual Position getTo() const = 0;

    /** return the length of the line */
    virtual Position getAlignedLength() const = 0;

    /** return the length of the line without gaps */
    virtual Position getFullLength() const = 0;

    /** add the specified number of gaps in the front and in the back */
    virtual void addGaps(int before, int after)	= 0;		
    
    /** add one or more gaps in the middle 
	@param position where gap(s) should be inserted 
	@param count the number of gaps to be inserted
    */
    virtual void insertGaps( int position, Position count = 1) = 0;
    
    /** remove leading/or trailing gaps */
    virtual void removeEndGaps() = 0;

    /** remove one or more positions from the aligned object */
    virtual void removeColumns( int position, Position count = 1) = 0; 

    /** remap the current alignatum object using ali 
	If unaligned chars is true, lower case unaligned characters will be
	put before the next aligned character (as much as fit) 
    */
    virtual void mapOnAlignment(const HAlignment & map_old2new,
				const Position new_length = 0,
				const bool unaligned_chars = false ) = 0;

    /** @short return an identical copy. */
    virtual HAlignatum getClone() const = 0;

    /** return an empty copy of this object */
    virtual HAlignatum getNew() const = 0;

    /** write into stream */
    virtual void  write( std::ostream & output ) const = 0;

    /** read from stream */
    virtual void  read( std::istream & input ) = 0;
    
};


}

#endif /* ALIGNATUM_H */

