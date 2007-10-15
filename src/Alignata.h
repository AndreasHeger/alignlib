/*
  alignlib - a library for aligning protein sequences

  $Id: Alignata.h,v 1.3 2004/03/19 18:23:39 aheger Exp $

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

#ifndef ALIGNATA_H
#define ALIGNATA_H 1

#include <iosfwd>
#include <string>
#include "alignlib.h"

namespace alignlib {

class Alignandum;
class AlignataIterator;
class AlignataConstIterator;

/**
   @short A residuepair containing row, column, and a score.
*/
struct ResiduePAIR {
    friend std::ostream & operator<< (std::ostream &, const ResiduePAIR &);

  ResiduePAIR() : mRow(0), mCol(0), mScore(0) {}
  ResiduePAIR( Position a, Position b, Score c) : mRow(a), mCol(b), mScore(c) {}
  ResiduePAIR& operator=( const ResiduePAIR & src) { mRow = src.mRow; mCol = src.mCol; mScore = src.mScore; return *this;}
  ResiduePAIR( const ResiduePAIR & src) : mRow(src.mRow), mCol(src.mCol), mScore(src.mScore) {}  
  Position mRow;
  Position mCol;
  Score mScore;
};

bool operator==( const ResiduePAIR & x, const ResiduePAIR & y); 
bool operator!=( const ResiduePAIR & x, const ResiduePAIR & y); 
 
/**     
    @short Protocol class for pairwise alignments.
	
    An alignment is a mapping.
    
    Each alignment should provide iterators. However,
    since the implementation differs between the different
    container types, I went for passive iterators, unlike
    the STL iterators, which are active. A passive iterator
    relies on the collection to manage its advancement and 
    therefore has to cooperate closely with its container.
    
    This class is a protocol class and as such defines only the
    general interface.
*/
 
class Alignata {
    friend std::ostream & operator<<(std::ostream &output, const Alignata &);
    friend std::istream & operator>>(std::istream &input, Alignata &);

 public:

    //------------------> constructors / destructors <---------------------------------------------------------
    /** empty constructor */
    Alignata();

    /** copy constructor */
    Alignata( const Alignata &src );

    /** destructor */
    virtual ~Alignata();

    //------------------------------------------------------------------------------------------------------------
    /** returns a new empty Alignment of the same type.
     */
    virtual Alignata * getNew() const = 0;
    
    /** returns an identical copy
     */
    virtual Alignata * getClone() const = 0;

    //------------------------------------------------------------------------------------------------------------
    /**
       @short Const iterator over an alignment.

       This iterator iterates over residues in the @ref Alignata
       objects.
     */
    class ConstIterator {
    public:
        ConstIterator() {};
 
        ConstIterator( const ConstIterator & src ) {};
	
        virtual ~ConstIterator() {};
	
	virtual ConstIterator * getClone() const = 0;

        /** dereference operator */
        virtual const ResiduePAIR & getReference() const = 0;

        /** for indirection */
        virtual const ResiduePAIR * getPointer() const = 0;
 
        /** advance one position */
        virtual void next() = 0;
	    
        /** step back one position */
	virtual void previous() = 0;
    };

    //------------------------------------------------------------------------------------------------------------
    /**
       @short Non-const iterator over an alignment.
     */
    class Iterator {
    public:
        Iterator() {};
 
        Iterator( const Iterator & src ) {};
 
        virtual ~Iterator() {};
	
	virtual Iterator * getClone() const = 0;

        /** dereference operator */
        virtual ResiduePAIR & getReference() const = 0;

        /** for indirection */
        virtual ResiduePAIR * getPointer() const = 0;
 
        /** advance one position */
        virtual void next() = 0;
	    
        /** step back one position */
	virtual void previous() = 0;
    };

    /** return const iterator */
    virtual AlignataConstIterator begin() const = 0; 
    
    /** return const iterator */
    virtual AlignataConstIterator end() const = 0; 

    /** return iterator */
    virtual AlignataIterator begin() = 0; 
    
    /** return iterator */
    virtual AlignataIterator end() = 0; 

    //----------------> accessors <------------------------------------------------------------------------------

    /** get the score of an alignment */
    virtual Score getScore() const = 0;

    /** get the length of an alignment */
    virtual Position	getLength() const = 0 ;

    /** get the number of gaps in the alignment */
    virtual Position getNumGaps() const = 0;

    /** set the score of an alignment */
    virtual void setScore( Score score ) = 0;

    /** returns true, if alignment is empty */
    virtual bool isEmpty() const = 0;

    /** returns the first residue aligned in row */
    virtual Position getRowFrom() const = 0;

    /** returns the last residue aligned in row */
    virtual Position	getRowTo() const = 0;

    /** returns the first residue aligned in col */
    virtual Position	getColFrom() const = 0;

    /** returns the last residue aligned in col */
    virtual Position	getColTo()   const = 0;

    /** returns the first aligned pair. Have to create a copy not a reference, because not all alignments will have
	a list of pairs ready */
    virtual ResiduePAIR getFirstPair() const = 0;
    
    /** returns the last aligned pair */
    virtual ResiduePAIR getLastPair() const = 0;

    /** adds a pair of residues to the alignment */
    virtual void addPair( ResiduePAIR * new_pair ) = 0; 

    /** adds a pair of residues to the alignment */
    virtual void addPair( Position row, Position col, Score score = 0) = 0; 

    /** removes a pair of residues from the alignment */
    virtual void removePair( const ResiduePAIR & old_pair ) = 0;

    /** retrieves a pair of residues from the alignment */
    virtual ResiduePAIR getPair( Position row ) const = 0;

    /** move alignment */
    virtual void moveAlignment( Position row_offset, Position col_offset) = 0;

    /** maps a residue from row to column. returns 0, if not found */
    virtual Position mapRowToCol( Position pos ) const = 0;

    /** maps a residue from row to column. returns 0, if not found */
    virtual Position mapColToRow( Position pos ) const = 0;

    /** removes a row-region in an alignment */
    virtual void removeRowRegion( Position from, Position to ) = 0;

    /** removes a column-region in an alignment */
    virtual void removeColRegion( Position from, Position to) = 0;

    /** switch row and column in the alignment */
    virtual void switchRowCol() = 0;

    /** clear the current alignemnt */
    virtual void clear() = 0;

    /*-----------------> I/O <------------------------------------------------------------------------------ */
      /** Write data members to stream */
    virtual void write(std::ostream & output ) const = 0;	       

    /** Read data members from stream */
    virtual void read( std::istream & input) = 0;
    
};

}
#endif /* _ALIGNATA_H */

