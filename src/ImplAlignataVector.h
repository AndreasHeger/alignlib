/*
  alignlib - a library for aligning protein sequences

  $Id: ImplAlignataVector.h,v 1.3 2004/03/19 18:23:40 aheger Exp $

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

#ifndef IMPL_ALIGNATA_VECTOR_H
#define IMPL_ALIGNATA_VECTOR_H 1

#include <iosfwd>
#include <vector>
#include "alignlib.h"
#include "ImplAlignata.h"

namespace alignlib {

class Alignandum;

/** alignment, where residues are kept in a vector.

    Residues are kept in a vector. 

    It is not assumed that the alignment is linear, meaning that
    for all i < j : col[i] < col[j] 

    The advantages:
    @li mapping is much faster

    The disadvantages:
    @li uses more memory for all the pairs, that are not part of the alignment
    @li iterating over all pairs is slower, since gaps have to be skipped explicitely.

    Implementation details:

    If the alignment is empty, mRowFrom is larger than mRowTo. This property is needed
    at some checks (for example in ImplAlignataVector_iterator)

    @author Andreas Heger
    @version $Id: ImplAlignataVector.h,v 1.3 2004/03/19 18:23:40 aheger Exp $
    @short protocol class for aligned objects
 */
class ImplAlignataVector : public ImplAlignata {

	typedef std::vector<ResiduePAIR *> PAIRVECTOR;
	typedef PAIRVECTOR::iterator PairIterator;
	typedef PAIRVECTOR::const_iterator PairConstIterator;

	friend class ConstIterator;

public:

	//------------------> constructors / destructors <---------------------------------------------------------
	/** constructor */
	ImplAlignataVector();

	/** copy constructor */
	ImplAlignataVector( const ImplAlignataVector &src );

	/** destructor */
	virtual ~ImplAlignataVector();

	//------------------------------------------------------------------------------------------------------------
	virtual ImplAlignataVector * getNew() const;

	/** return an identical copy */
	virtual ImplAlignataVector * getClone() const;

	//------------------------------------------------------------------------------------------------------------
	//------------------------------------------------------------------------------------------------------------
	/**
       @short Const iterator over an alignment.
	 */
	class ImplAlignataVector_ConstIterator : public Alignata::ConstIterator{
	public:
		ImplAlignataVector_ConstIterator(const PAIRVECTOR & container, 
				Position current, 
				Position from, 
				Position to) : 
					mContainer( container ),
					mCurrentRow( current ), 
					mFirstRow( from), 
					mLastRow( to ) { 
			if (mCurrentRow >= mLastRow || to == NO_POS) 
				mCurrentRow = NO_POS;
		};

		ImplAlignataVector_ConstIterator( const ImplAlignataVector_ConstIterator & src ) : 
			mContainer( src.mContainer),
			mCurrentRow( src.mCurrentRow), 
			mFirstRow( src.mFirstRow), 
			mLastRow( src.mLastRow) {};

			virtual ~ImplAlignataVector_ConstIterator() {};

			virtual ConstIterator * getClone() const {return new ImplAlignataVector_ConstIterator( *this );}

			/** dereference operator: runtime error, if out of bounds? */
			virtual const ResiduePAIR & getReference() const { return *mContainer[mCurrentRow];}

			/** for indirection */
			virtual const ResiduePAIR * getPointer() const { 
				if (mCurrentRow != NO_POS) 
					return mContainer[mCurrentRow]; 
				else 
					return NULL;
			}

			/** advance one position, until you find an aligned pair */
			virtual void next() { 
				mCurrentRow++;
				while (mCurrentRow < mLastRow && mContainer[mCurrentRow] == NULL) { 
					mCurrentRow++;
				}
				if (mCurrentRow >= mLastRow) mCurrentRow = NO_POS;
			}

			/** step back one position, until you find an aligned pair */
			virtual void previous() { 
				mCurrentRow--; 
				while (mCurrentRow >= mFirstRow && mContainer[mCurrentRow] == NULL) { 
					mCurrentRow--;
				}
				if (mCurrentRow < mFirstRow) mCurrentRow = NO_POS;
			}

	private:
		const PAIRVECTOR & mContainer;
		Position mCurrentRow;
		Position mFirstRow;
		Position mLastRow;
	};

	//------------------------------------------------------------------------------------------------------------
	/**
       @short Non-Const iterator over an alignment.
	 */
	class ImplAlignataVector_Iterator : public Alignata::Iterator{
	public:
		ImplAlignataVector_Iterator(PAIRVECTOR & container,
				Position current,
				Position from,
				Position to) : 
					mContainer( container ),
					mCurrentRow( current ), 
					mFirstRow( from), 
					mLastRow( to ) { 
			if (mCurrentRow >= mLastRow || to == NO_POS) 
				mCurrentRow = NO_POS;
		};

		ImplAlignataVector_Iterator( const ImplAlignataVector_Iterator & src ) : 
			mContainer( src.mContainer),
			mCurrentRow( src.mCurrentRow), 
			mFirstRow( src.mFirstRow), 
			mLastRow( src.mLastRow) {};

			virtual ~ImplAlignataVector_Iterator() {};

			virtual Iterator * getClone() const {return new ImplAlignataVector_Iterator( *this );}

			/** dereference operator: runtime error, if out of bounds? */
			virtual ResiduePAIR & getReference() const { return *mContainer[mCurrentRow];}

			/** for indirection */
			virtual ResiduePAIR * getPointer() const { 
				if (mCurrentRow != NO_POS) 
					return mContainer[mCurrentRow]; 
				else 
					return NULL;
			}

			/** advance one position, until you find an aligned pair */
			virtual void next() { 
				mCurrentRow++;
				while (mCurrentRow < mLastRow && mContainer[mCurrentRow] == NULL ) { 
					mCurrentRow++;
				}
				if (mCurrentRow >= mLastRow) mCurrentRow = NO_POS;
			}

			/** step back one position, until you find an aligned pair */
			virtual void previous() { 
				mCurrentRow--; 
				while (mCurrentRow >= mFirstRow && mContainer[mCurrentRow] == NULL) { 
					mCurrentRow--;
				}
				if (mCurrentRow < mFirstRow) mCurrentRow = NO_POS;
			}

	private:
		PAIRVECTOR & mContainer;
		Position mCurrentRow;
		Position mFirstRow;
		Position mLastRow;
	};

	/** return const iterator */
	virtual AlignataConstIterator begin() const; 

	/** return const iterator */
	virtual AlignataConstIterator end() const; 

	/** return const iterator */
	virtual AlignataIterator begin(); 

	/** return const iterator */
	virtual AlignataIterator end(); 

	//----------------> accessors <------------------------------------------------------------------------------

	/** returns the first residue aligned in row */
	virtual Position getRowFrom() const;

	/** returns the last residue aligned in row */
	virtual Position	getRowTo() const;

	/** returns the first residue aligned in col */
	virtual Position	getColFrom() const;

	/** returns the last residue aligned in col */
	virtual Position	getColTo()   const;

	/** returns the first aligned pair. Have to create a copy not a reference, because not all alignments will have
	a list of pairs ready */
	virtual ResiduePAIR getFirstPair() const;

	/** returns the last aligned pair */
	virtual ResiduePAIR getLastPair() const;

	/** adds a pair of residue to the alignment */
	virtual void addPair( ResiduePAIR * new_pair ); 

	/** removes a pair of residues from the alignment */
	virtual void removePair( const ResiduePAIR & old_pair );

	/** retrieves a pair of residues from the alignment */
	virtual ResiduePAIR getPair( Position row ) const;

	/** maps a residue from row to column. returns 0, if not found. This is quick, since there is 
	only one lookup in the array needed.*/
	virtual Position mapRowToCol( Position pos, SearchType search = NO_SEARCH ) const;

	/** move alignment */
	virtual void moveAlignment( Position row_offset, Position col_offset);

	/** removes a row-region in an alignment */
	virtual void removeRowRegion( Position from, Position to );

	/** removes a column-region in an alignment */
	virtual void removeColRegion( Position from, Position to);

	/** clear the current alignemnt */
	void clear();
	protected:

		/** reset boundaries mRowFrom and mRowTo after the deletion of pairs */
		void resetBoundaries();

	private:

		/** clear container */
		void clearContainer();

		/** Vector of residue pairs */
		PAIRVECTOR mPairs;

		/** row of the first aligned residue pair */
		Position mRowFrom;

		/** row of the last aligned residue pair */
		Position mRowTo;
};



}

#endif /* _ALIGNATA_H */

