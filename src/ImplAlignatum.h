/*
  alignlib - a library for aligning protein sequences

  $Id: ImplAlignatum.h,v 1.5 2004/03/19 18:23:41 aheger Exp $

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

#ifndef IMPL_ALIGNATUM_H
#define IMPL_ALIGNATUM_H 1

#include <iosfwd>
#include <string>
#include "alignlib.h"
#include "alignlib_fwd.h"
#include "Alignatum.h"

namespace alignlib 
{

/** @short Implementation class of @ref Alignatum objects.

    @author Andreas Heger
    @version $Id: ImplAlignatum.h,v 1.5 2004/03/19 18:23:41 aheger Exp $
*/
class ImplAlignatum : public Alignatum 
{
 public:
    // constructors and desctructors
    /** when building a new Alignatum sequence by cloning and an alignment, this specifies,
	if the new sequence is in row or column
    */

    /** constructor */
    ImplAlignatum  ();
    
    /** copy constructor */
    ImplAlignatum  (const ImplAlignatum & src);
    /** destructor */
    virtual ~ImplAlignatum ();

    ImplAlignatum (const std::string & representation, 
    		Position from = NO_POS, 
    		Position to = NO_POS);

    /*-----> accessors <----------------------------------------------------- */
    /** write object into stream nicely formatted. Segments are addressed by [from,to)
	@param output stream for result
	@param segment_start beginning of segment 
	@param segment_end end of segment 
    */
    virtual void writeRow( std::ostream & output,
    		const HRenderer & renderer,
			Position segment_start = NO_POS, 
			Position segment_end = NO_POS ) const;

    /** readline */
    virtual void readRow( std::istream & input );

    /** return a copy of the string representation */
    virtual std::string getString() const;

    /** return a reference to the string representation */
    virtual const std::string & getStringReference() const;

    /** return the length of the line */
    virtual Position getAlignedLength() const;

    /** return the length of the line without gaps */
    virtual Position getFullLength() const;
    
    /** return the number of the first residue */
    virtual Position getFrom() const;

    /** return the number of the last residue */
    virtual Position getTo() const;

    /** add the specified number of gaps in the front and in the back */
    virtual void addGaps(int before, int after);		
    
    /** add one or more gaps in the middle 
	@param position where gap(s) should be inserted 
	@param count the number of gaps to be inserted
    */
    virtual void insertGaps( int position, Position count = 1);
    
    /** remove leading/or trailing gaps */
    virtual void removeEndGaps();

    /** remove one or more positions from the aligned object */
    virtual void removeColumns( int position, Position count = 1); 

    /** remap the current alignatum object using ali. 
	If unaligned chars is true, lower case unaligned characters will be
	put before the next aligned character (as much as fit) 
    */
    virtual void mapOnAlignment(
    		const HAlignment & map_old2new,
    		const Position new_length = 0,
    		const bool unaligned_chars = false );

    /** return a copy of the object, mapped if so desired 
    */
    virtual HAlignatum getClone() const;
    
    /** return an empty copy of this object */
    virtual HAlignatum getNew() const;

    /** write into stream */
    virtual void  write( std::ostream & output ) const;

    /** read from stream */
    virtual void  read( std::istream & input );

 protected:

    /** calculate number of gaps in sequence */
    virtual int countGaps();

    /** calculate residue number of residue in position pos, or, if this is a gap, the residue number
	of the next non-gap character. */
    virtual Position getResidueNumberNext( Position pos ) const;

    /** calculate residue number of residue in pos, or, if this is a gap, the residue number
	of the previous non-gap character. */
    virtual Position getResidueNumberPrevious( Position pos ) const;

    /** accessor functions */
    
    /** get represenation */
    const std::string & getRepresentation() const;

    /** set representation */
    void setRepresentation( std::string & representation,
			    Position first_res = NO_POS,
			    Position last_res = NO_POS);

    inline char getFieldSeparator() const { return mSeparator; }

 private:
    
    /** the string representation used for single line output */
    std::string mRepresentation;

    /** number of first residue in aligned string */
    Position mFrom;
    
    /** number of last residue in aligned string */
    Position mTo;
    
    /** length of aligned string */
    Position mLength;

    /** gap-character */
    char mGapChar;

    /** field separator */
    char mSeparator;
};

}

#endif /* ALIGNATUM_H */

