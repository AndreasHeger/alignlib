/*
  alignlib - a library for aligning protein sequences

  $Id: ImplAlignataHash.h,v 1.3 2004/03/19 18:23:40 aheger Exp $

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

#ifndef IMPL_ALIGNATA_HASH_H
#define IMPL_ALIGNATA_HASH_H 1

#include <iosfwd>
#include "alignlib.h"
#include "ImplAlignata.h"

namespace alignlib
{

class Alignandum;

/**
    @short alignment, where residues are stored in a hash sorted by row and then col
    
    Residues are kept in a set, but in contrast to AlignataSet, both row and
    column are taken into account, when comparing keys.
    
    The code is basically identical to ImplAlignataSet, the comparator has 
    changed. Maybe in the future this can be made nicer by subclassing or
    template instantiation.
	
    @author Andreas Heger
    @version $Id: ImplAlignataHash.h,v 1.3 2004/03/19 18:23:40 aheger Exp $
*/
class ImplAlignataHash : public ImplAlignata {

  /** residues are sorted by row and then by column */
  struct Comparator {
      bool operator() ( const ResiduePAIR * x, const ResiduePAIR * y) const { 
	  if (x->mRow < y->mRow) return 1;
	  if (x->mRow > y->mRow) return 0;
	  if (x->mCol < y->mCol) 
	      return 1;
	  else
	      return 0;
      }
  };  
  typedef std::set<ResiduePAIR *, Comparator> PAIRSET;
  typedef PAIRSET::const_iterator PairConstIterator;
  typedef PAIRSET::iterator PairIterator;
  
 public:
    //------------------> constructors / destructors <---------------------------------------------------------
    /** constructor */
    ImplAlignataHash();

    /** copy constructor */
    ImplAlignataHash( const ImplAlignataHash &src );

    /** destructor */
    virtual ~ImplAlignataHash();

    //------------------------------------------------------------------------------------------------------------
    virtual ImplAlignataHash * getNew() const;
    
    /** return an identical copy */
    virtual ImplAlignataHash * getClone() const;

    //------------------------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------------------------
    
    /**
       @short Const iterator over an alignment.

       Const iterator over alignment residue pairs.
     */
    
    class ImplAlignataHash_ConstIterator : public Alignata::ConstIterator {
    public:
        ImplAlignataHash_ConstIterator(PairConstIterator it) : mIterator( it ) {};
	
        ImplAlignataHash_ConstIterator( const ImplAlignataHash_ConstIterator & src ) : mIterator( src.mIterator) {};
	
        virtual ~ImplAlignataHash_ConstIterator() {};
	
  	virtual ImplAlignataHash_ConstIterator * getClone() const {return new ImplAlignataHash_ConstIterator( mIterator);}

        /** dereference operator */
        virtual const ResiduePAIR & getReference() const { return *(*mIterator);}

	/** for indirection */
        virtual const ResiduePAIR * getPointer() const { return (*mIterator); }

        /** advance one position */
        virtual void next() { mIterator++;}
	
        /** step back one position */
	virtual void previous() { mIterator--;}
    private:
	PairConstIterator mIterator;
    };

    //------------------------------------------------------------------------------------------------------------
    /**
       @short Non-Const iterator over an alignment.
     */
    class ImplAlignataHash_Iterator : public Alignata::Iterator {
    public:
        ImplAlignataHash_Iterator(PairIterator it) : mIterator( it ) {};
	
        ImplAlignataHash_Iterator( const ImplAlignataHash_Iterator & src ) : mIterator( src.mIterator) {};
	
        virtual ~ImplAlignataHash_Iterator() {};
	
  	virtual ImplAlignataHash_Iterator * getClone() const {return new ImplAlignataHash_Iterator( mIterator);}

        /** dereference operator
	    The solution here compiles at least. It is ugly and inefficient. However, when returning *mIterator, you get the following error:
	    "alignlib::ResiduePAIR &" to an initializer of type "const
	    std::_Rb_tree<std::set<alignlib::ResiduePAIR,               
	    ...
	 */
        virtual ResiduePAIR & getReference() const { return *(*mIterator);}

	/** for indirection 
	    The solution here compiles at least. When returning return (mIterator.operator->()), the const-type is not correctly
	    set.
	 */
        virtual ResiduePAIR * getPointer() const { return *mIterator; }	

        /** advance one position */
        virtual void next() { mIterator++; }
	    
        /** step back one position */
	virtual void previous() { mIterator--;}
    private:
	PairIterator mIterator;
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

    /** retrieves a pair of residues from the alignment */
    virtual ResiduePAIR getPair( Position row ) const;

    /** removes a pair of residues from the alignment */
    virtual void removePair( const ResiduePAIR & old_pair );

    /** clear the current alignemnt */
    void clear();

 private:
     PAIRSET mPairs;
};

						  

}

#endif /* _ALIGNATA_H */

