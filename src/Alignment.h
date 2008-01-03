/*
  alignlib - a library for aligning protein sequences

  $Id: Alignment.h,v 1.3 2004/03/19 18:23:39 aheger Exp $

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

namespace alignlib 
{

class Alignandum;
class AlignmentIterator;
class AlignmentConstIterator;

/**
   @short A residuepair containing row, column, and a score.
 */
struct ResiduePAIR 
{
	/** Base class for residue aligned residues.
	 */

	friend std::ostream & operator<< (std::ostream &, const ResiduePAIR &);

	ResiduePAIR() : mRow(NO_POS), mCol(NO_POS), mScore(0) {}
	ResiduePAIR( Position a, Position b, Score c = 0) : mRow(a), mCol(b), mScore(c) {}
	ResiduePAIR& operator=( const ResiduePAIR & src) 
	{ 
		mRow = src.mRow; 
		mCol = src.mCol; 
		mScore = src.mScore; 
		return *this;
	}
	ResiduePAIR( const ResiduePAIR & src) : mRow(src.mRow), mCol(src.mCol), mScore(src.mScore) {}  
	Position mRow;
	Position mCol;
	Score mScore;
};

bool operator==( const ResiduePAIR & x, const ResiduePAIR & y); 
bool operator!=( const ResiduePAIR & x, const ResiduePAIR & y); 

/**     
    @short Protocol class for pairwise alignments.

    Tthe purpose of this class to provide a mapping between
    two objects. In true dynamic programming spirit, the
    two objects are called "row" and "col".

    Different implementations of this interface store residues
    in different order and implement different efficiences in
    mapping residues.

    Each alignment provides iterators. However,
    since the implementation differs between the different
    container types, I went for passive iterators, unlike
    the STL iterators, which are active. A passive iterator
    relies on the collection to manage its advancement and 
    therefore has to cooperate closely with its container.

    This class is a protocol class and as such defines only 
    the general interface.


 */

class Alignment 
{
	friend std::ostream & operator<<(std::ostream &output, const Alignment &);

public:

	//------------------> constructors / destructors <---------------------------------------------------------
	/** empty constructor */
	Alignment();

	/** copy constructor */
	Alignment( const Alignment &src );

	/** destructor */
	virtual ~Alignment();

	//------------------------------------------------------------------------------------------------------------
	/** returns a new empty Alignment of the same type.
	 */
	virtual Alignment * getNew() const = 0;

	/** returns an identical copy
	 */
	virtual Alignment * getClone() const = 0;

	//------------------------------------------------------------------------------------------------------------
	/**
       @short Const iterator over an alignment.

       This iterator iterates over residues in the @ref Alignment
       objects.
	 */
	class ConstIterator 
	{
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
	class Iterator 
	{
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

	/** returns a const_iterator pointing to the first element in the container
	 */
	virtual AlignmentConstIterator begin() const = 0; 

	/** returns a const_iterator pointing one past the last element in the container.
	 */
	virtual AlignmentConstIterator end() const = 0; 

	/** returns an iterator pointing to the first element in the container
	 */
	virtual AlignmentIterator begin() = 0; 

	/** returns an iterator pointing one past the last element in the container.
	 */
	virtual AlignmentIterator end() = 0; 

	//----------------> accessors <------------------------------------------------------------------------------

	/** get the score of an alignment
	 */
	virtual Score getScore() const = 0;

	/** get the length of an alignment 

     The length of an alignment is the number of aligned residue pairs plus
     the number gaps.
	 */
	virtual Position getLength() const = 0 ;

	/** get the number of gaps in the alignment 
	 * 
	 * Gaps are counted in either sequence.
	 * */
	virtual Position getNumGaps() const = 0;

    /** get the number of aligned positions 
     * */
    virtual Position getNumAligned() const = 0;
	
	/** set the alignment score. 
	 * */
	virtual void setScore( Score score ) = 0;

	/** returns true, if alignment is empty */
	virtual bool isEmpty() const = 0;

	/** returns the first residue aligned in row 
	 * */
	virtual Position getRowFrom() const = 0;

	/** returns one past the last residue aligned in row 
	 * */
	virtual Position getRowTo() const = 0;

	/** returns the first residue aligned in col 
	 * */
	virtual Position getColFrom() const = 0;

	/** returns one past the last residue aligned in col 
	 * */
	virtual Position getColTo() const = 0;

	/** returns a copy of the first aligned pair. 
	 * */
	virtual ResiduePAIR front() const = 0;

	/** returns a copy of the last aligned pair 
	 * */
	virtual ResiduePAIR back() const = 0;

	/** adds a pair of residues to the alignment 
	 * 
	 * The alignata object takes ownership of the residue pair.
	 * */
	virtual void addPair( ResiduePAIR * new_pair ) = 0; 

	/** adds a pair of residues to the alignment */
	virtual void addPair( Position row, Position col, Score score = 0) = 0; 

	/** removes a pair of residues from the alignment */
	virtual void removePair( const ResiduePAIR & old_pair ) = 0;

	/** retrieves a copy of a residue pair in row from the 
	 * alignment
	 *  
	 * @param row row of residue pair
	 * */
	virtual ResiduePAIR getPair( const ResiduePAIR & p) const = 0;

	/** move alignment */
	virtual void moveAlignment( Position row_offset, Position col_offset) = 0;

	/** maps a residue from row to column. returns 0, if not found */
	virtual Position mapRowToCol( Position pos, SearchType search = NO_SEARCH ) const = 0;

	/** maps a residue from row to column. returns 0, if not found */
	virtual Position mapColToRow( Position pos, SearchType search = NO_SEARCH ) const = 0;

	/** removes a row-region in an alignment */
	virtual void removeRowRegion( Position from, Position to ) = 0;

	/** removes a column-region in an alignment */
	virtual void removeColRegion( Position from, Position to) = 0;

	/** switch row and column in the alignment */
	virtual void switchRowCol() = 0;

	/** clear the current alignemnt */
	virtual void clear() = 0;

	/*-----------------> I/O <------------------------------------------------------------------------------ */
	/** write human readable representation of alignment to stream
	 */
	virtual void write(std::ostream & output ) const = 0;	       

};

}
#endif /* _ALIGNATA_H */

