/*
  alignlib - a library for aligning protein sequences

  $Id: ImplAlignata.h,v 1.3 2004/03/19 18:23:40 aheger Exp $

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

#ifndef IMPL_ALIGNATA_H
#define IMPL_ALIGNATA_H 1

#include <iosfwd>
#include "alignlib.h"
#include "Alignata.h"

namespace alignlib
{

  /** basic implementation for aligned objects. 

      @author Andreas Heger
      @version $Id: ImplAlignata.h,v 1.3 2004/03/19 18:23:40 aheger Exp $
      @short basic implementation class for aligned objects
  */

class ImplAlignata : public Alignata {

 public:
    //------------------> constructors / destructors <---------------------------------------------------------
    /** empty constructor */
    ImplAlignata();

    /** copy constructor */
    ImplAlignata( const ImplAlignata &src );

    /** destructor */
    virtual ~ImplAlignata();
    
    //----------------> accessors <------------------------------------------------------------------------------
    virtual Score getScore() const;

    virtual Position	getLength() const;

    virtual Position getNumGaps() const;

    virtual void setScore( Score score );

    bool isEmpty() const;

    virtual void switchRowCol();

    virtual void moveAlignment( Position row_offset, Position col_offset);

    /** maps a residue from row to column. returns 0, if not found. This default, but working implementation is
	very time-inefficient, especially, if you want to map the whole thing. Other implementations in derived
	classes can be much faster */
    virtual Position mapRowToCol( Position pos ) const;

    /** maps a residue from column to row. returns 0, if not found. This default, but working implementation is
	very time-inefficient, especially, if you want to map the whole thing. Other implementations in derived
	classes can be much faster*/
    virtual Position mapColToRow( Position pos ) const;

    virtual void removeRowRegion( Position from, Position to );

    virtual void removeColRegion( Position from, Position to);

    virtual void clear() ;

    /** @brief adds a pair of residues to the alignment 
	(have to add this here, otherwise it won't compile!. It seems that overloaded
	functions can not be separately implemented) */
    virtual void addPair( ResiduePAIR * new_pair ) = 0; 

    /** adds a pair of residues to the alignment */
    virtual void addPair( Position row, Position col, Score score = 0); 
    
    /*-----------------> I/O <------------------------------------------------------------------------------ */
    virtual void write(std::ostream & output ) const;	       

    virtual void read( std::istream & input);

 protected:
    /** the length has changed of the alignment */
    virtual void setChangedLength();

    /** set length of the alignment */
    virtual void setLength( Position length ) const;

    /** set length of the alignment */
    virtual void setNumGaps( Position num_gaps ) const;

    /** calculate alignment length */
    virtual void calculateLength() const;

    /** flag, whether alignment length has changed. Protected, because some childs need to query and 
     set this flag.
    */
    mutable bool mChangedLength;

 private:

    /** length of alignemnt, has to be mutable, since length evaluation is lazy */
    mutable Position mLength;
    
    /** total score of alignment */
    Score mScore;
    
    /** number of gaps in alignment */
    mutable Position mNumGaps;

    
};
						  

}

#endif /* IMPL_ALIGNATA_H */

