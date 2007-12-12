/*
  alignlib - a library for aligning protein sequences

  $Id: ImplAlignataSetCol.h,v 1.3 2004/03/19 18:23:40 aheger Exp $

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

#ifndef IMPL_ALIGNATA_SET_COL_H
#define IMPL_ALIGNATA_SET_COL_H 1

#include <iosfwd>
#include "alignlib.h"
#include "ImplAlignata.h"

namespace alignlib {

class Alignandum;
 
/** @short alignment, where residues are stored in a hash sorted by col and then by row.
	
    @author Andreas Heger
    @version $Id: ImplAlignataSetCol.h,v 1.3 2004/03/19 18:23:40 aheger Exp $
    @short protocol class for aligned objects
*/
class ImplAlignataSetCol : public ImplAlignata {

  /** residues are sorted by row */
  struct Comparator {
    bool operator() ( const ResiduePAIR * x, const ResiduePAIR * y) const { return x->mCol < y->mCol; } 
  };
    
  typedef std::set<ResiduePAIR *, Comparator> PAIRSET;
  typedef PAIRSET::const_iterator PairConstIterator;
  typedef PAIRSET::iterator PairIterator;
  
 public:

    //------------------> constructors / destructors <---------------------------------------------------------
    /** constructor */
    ImplAlignataSetCol();

    /** copy constructor */
    ImplAlignataSetCol( const ImplAlignataSetCol &src );

    /** destructor */
    virtual ~ImplAlignataSetCol();

    //------------------------------------------------------------------------------------------------------------
    virtual ImplAlignataSetCol * getNew() const;
    
    /** return an identical copy */
    virtual ImplAlignataSetCol * getClone() const;

    //------------------------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------------------------
    /**
       @short Const iterator over an alignment.
     */
    class ImplAlignataSetCol_ConstIterator : public Alignata::ConstIterator {
    public:
        ImplAlignataSetCol_ConstIterator(PairConstIterator it) : mIterator( it ) {};
	
        ImplAlignataSetCol_ConstIterator( const ImplAlignataSetCol_ConstIterator & src ) : mIterator( src.mIterator) {};
	
        virtual ~ImplAlignataSetCol_ConstIterator() {};
	
  	virtual ImplAlignataSetCol_ConstIterator * getClone() const {return new ImplAlignataSetCol_ConstIterator( mIterator);}

        /** dereference operator */
        virtual const ResiduePAIR & getReference() const { return **mIterator;}

	/** for indirection */
        virtual const ResiduePAIR * getPointer() const { return *mIterator; }

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
    class ImplAlignataSetCol_Iterator : public Alignata::Iterator {
    public:
        ImplAlignataSetCol_Iterator(PairIterator it) : mIterator( it ) {};
	
        ImplAlignataSetCol_Iterator( const ImplAlignataSetCol_Iterator & src ) : mIterator( src.mIterator) {};
	
        virtual ~ImplAlignataSetCol_Iterator() {};
	
  	virtual ImplAlignataSetCol_Iterator * getClone() const {return new ImplAlignataSetCol_Iterator( mIterator);}

        /** dereference operator
	    The solution here compiles at least. It is ugly and inefficient. However, when returning *mIterator, you get the following error:
	    "alignlib::ResiduePAIR &" to an initializer of type "const
	    std::_Rb_tree<std::set<alignlib::ResiduePAIR,               
	    ...
	 */
        virtual ResiduePAIR & getReference() const { return **mIterator;}

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
    virtual ResiduePAIR front() const;
    
    /** returns the last aligned pair */
    virtual ResiduePAIR back() const;

    /** adds a pair of residue to the alignment */
    virtual void addPair( ResiduePAIR * new_pair ); 

    /** removes a pair of residues from the alignment */
    virtual void removePair( const ResiduePAIR & old_pair );

    /** retrieves a pair of residues from the alignment */
    virtual ResiduePAIR getPair( const ResiduePAIR & p) const;

    /** clear the current alignemnt */
    virtual void clear();

    /** removes a row-region in an alignment */
    virtual void removeRowRegion( Position from, Position to );

    /** removes a row-region in an alignment */
    virtual void removeColRegion( Position from, Position to );
    
 protected:

 private:

    /** clear container */
    void clearContainer();
    
    /** container with residue pairs */
    PAIRSET mPairs;
};

						  

}

#endif /* _ALIGNATA_H */

