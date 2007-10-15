/*
  alignlib - a library for aligning protein sequences

  $Id: ImplAlignandum.h,v 1.2 2004/01/07 14:35:33 aheger Exp $

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

#ifndef IMPL_ALIGNANDUM_H
#define IMPL_ALIGNANDUM_H 1

#include <iosfwd>
#include "Alignandum.h"

namespace alignlib {

/** 
    Base class for objects that are to be aligned. This class implements a subset of the interface,
    namely some basic behaviour for using segments and communicating its state (prepared/not prepared), 
    that is common to all its derived classes.

    Translator: There is one translator-object in the whole class library. If an object needs a different
    translator, then the default translator has to be changed. See Translator for more information.

    @author Andreas Heger
    @version $Id: ImplAlignandum.h,v 1.2 2004/01/07 14:35:33 aheger Exp $
    @short protocol class of alignable objects
*/

class ImplAlignandum : public Alignandum {
    /* friends ---------------------------------------------------------------------------- */
    friend  std::ostream & operator<<( std::ostream &, const ImplAlignandum &);
    friend  std::istream & operator>>( std::istream &, ImplAlignandum &);             
    
 public:
    /* constructors and desctructors------------------------------------------------------- */
  
    /** empty constructor */
    ImplAlignandum();
    
    /** copy constructor */
    ImplAlignandum( const ImplAlignandum &);

    /** desctructor */
    virtual ~ImplAlignandum();

    /** accessors ------------------------------------------------------------------------- */
    /** get length of window */
    virtual Position	getLength() const;
    
    /** use full length of sequence for exporting. Basically extends the region to the
	true begining and end of object */
    virtual void useFullLength();
    
    /** use a segment for exporting and set segment to from and to 
	@param from	where segment starts
	@param to		where segment ends
    */
    virtual void useSegment( Position from, Position to);

    /** return true if object is prepared for alignment (for cacheable types ) */
    virtual bool isPrepared() const;	

    /** get internal representation of residue in position pos */
    virtual char asChar( Position pos ) const;

    /** returns a string representation of the object */
    virtual std::string asString() const;

    /** mask column in a segment */
    virtual void mask( Position from, Position to);
    
    /** mask column in a segment */
    virtual void mask( Position column) = 0;

    /** return first residue number in segment */
    virtual Position getFrom() const;

    /** return last residue number in segment */
    virtual Position getTo() const;

 protected:
    /** the member functions below are protected, because they have to be only accessible for
	derived classes. They should know, what they are doing. */


    /** given a position x in a segment, return the offset from the true beginning of the sequence */
    virtual Position getOffset( Position x ) const;
    
    /** set true length*/
    virtual void setTrueLength(Position length) const;

    /** get true length*/
    virtual Position getTrueLength() const;

    /** set prepared flag */
    virtual void setPrepared( bool flag ) const;

 private:
    /** first residue of segment used for aligning */
    Position  mFrom;

    /** last residue of segment used for aligning */
    Position  mTo;                   

    /** true length of sequence */
    mutable Position mLength; 

    /** flag, whether object is ready for alignment */
    mutable bool mIsPrepared;                          

};


}

#endif /* ALIGNANDUM_H */

