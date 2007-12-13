/*
  alignlib - a library for aligning protein sequences

  $Id: ImplSequence.h,v 1.2 2004/01/07 14:35:36 aheger Exp $

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

#ifndef IMPL_SEQUENCE_H
#define IMPL_SEQUENCE_H 1

#include "alignlib.h"
#include "ImplAlignandum.h"
#include <iosfwd>


namespace alignlib {

    class Alignment;

struct AlignandumDataSequence : public AlignandumData {
  Residue * mSequencePointer;
};
    /** A class for sequences, that are to be aligned. Instances of this
	class are created by factory functions. This class implements the
	part of the interface, that has not been implemented by IAlignandum
	
	@author Andreas Heger
	@version $Id: ImplSequence.h,v 1.2 2004/01/07 14:35:36 aheger Exp $
	@short Contains sequence, that are to be aligned
     */


class ImplSequence : public ImplAlignandum {

    friend Alignandum * addSequence2Profile( Alignandum * dest, const Alignandum * source, const Alignment * map_source2dest );

 public:
    /*------------------------------------------------------------------------------------ */
    /** Since I use lazy evaluation/retrieval of sequences from databases, I need an 
	empty constructor */
    ImplSequence();

    /** create sequence from a string, given a translator object */
    ImplSequence( const char * src );

    /** the copy constructor */
    ImplSequence( const ImplSequence & );

    /** the destructor */
    virtual ~ImplSequence();
    
    /*------------------------------------------------------------------------------------ */
    /** return an identical copy of this object */
    virtual Alignandum * getClone() const;

    /** return a pointer to the member data of this sequence */
    virtual const AlignandumDataSequence & getData() const;

    /** get internal representation of residue in position pos */
    virtual Residue asResidue( Position pos ) const;
    
    /** mask column at position x */
    virtual void mask( Position x);

    /** shuffle object */
    virtual void shuffle( unsigned int num_iterations = 1,
			  Position window_size = 0);

    /* Mutators ------------------------------------------------------------------------------ */
    
    /** load data into cache, if cacheable type */
    virtual void prepare() const;						
    
    /** discard cache, if cacheable type */
    virtual void release() const;					       

    /** write human readable output to stream.
     */
    virtual void write( std::ostream & output ) const;
    
    /** save state of object into stream
     */
    virtual void load( std::istream & input ) ;
    
    
 protected:

	 /** save state of object into stream
	  */
	 virtual void __save( std::ostream & output, MagicNumberType type = MNNoType ) const;

 private:

 public:

 protected:
    /** where I store my member data */
    mutable AlignandumDataSequence mData;				

    /** The actual sequence is here*/
    Residue * mSequence;
};

}

#endif /* _SEQUENCE_H */

