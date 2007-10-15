/*
  alignlib - a library for aligning protein sequences

  $Id: AlignataIterator.h,v 1.3 2004/03/19 18:23:39 aheger Exp $

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

#ifndef ALIGNATA_ITERATOR_H
#define ALIGNATA_ITERATOR_H 1

#include <iosfwd>
#include "alignlib.h"
#include "Alignata.h"

namespace alignlib {

    /** This class takes any of the iterators
	of aligned pairs of alignata-objects
	and converts it into an iterator, that
	is STL-compatible in its syntax (i.e. 
	no overloading, etc.). The penalty for
	this universality in syntax is paid in
	performance. There are two additional
	indirections.

	@author Andreas Heger
	@version $Id: AlignataIterator.h,v 1.3 2004/03/19 18:23:39 aheger Exp $
	@short protocol class for aligned objects
    */

class AlignataConstIterator {
 public:

    AlignataConstIterator() : mIterator(NULL) {};

    AlignataConstIterator( Alignata::ConstIterator * it) : mIterator(it) {} ;
    
    AlignataConstIterator( const AlignataConstIterator & src ) { mIterator = src.mIterator->getClone();};
    
    ~AlignataConstIterator() { delete mIterator; };

    /** assignment operator */
    inline AlignataConstIterator & operator=( const AlignataConstIterator & other ) {
      if (*this != other) {
	/* mIterator = other.mIterator; does not work, I guess, because of shallow copying */
	delete mIterator;
	mIterator = other.mIterator->getClone();
      }
      
      return *this;
    }

    /** comparison operator */
    inline bool operator==( const AlignataConstIterator & other) {   
      return ( (*mIterator).getPointer() == (*other.mIterator).getPointer() );
    }
    
    /** comparison operator */
    inline bool operator!=( const AlignataConstIterator & other) {
      return ( (*mIterator).getPointer() != (*other.mIterator).getPointer() );
    }
      
    /** indirection operator */
    inline const ResiduePAIR * operator->() const { return &(mIterator->getReference()); }

    /** dereference operator */
    inline const ResiduePAIR & operator*() const { return mIterator->getReference(); }
 
    /** prefix ++ */
    inline AlignataConstIterator & operator++() { mIterator->next(); return *this; }
 
    /** postfix ++ */
    inline AlignataConstIterator operator++(int) { AlignataConstIterator tmp = *this; mIterator->next(); return tmp; }
 
    /** prefix -- */
    inline AlignataConstIterator & operator--() { mIterator->previous(); return *this; }
 
    /** postfix -- */
    inline AlignataConstIterator  operator--(int) { 
      AlignataConstIterator tmp = *this; mIterator->previous(); return tmp; }
    
 private:
    Alignata::ConstIterator * mIterator;
};

class AlignataIterator {
 public:
    AlignataIterator() : mIterator(NULL) {};

    AlignataIterator( Alignata::Iterator * it) : mIterator(it) {} ;
    
    AlignataIterator( const AlignataIterator & src ) { mIterator = src.mIterator->getClone();};
    
    ~AlignataIterator() { delete mIterator; };

    /** assignment operator */
    inline AlignataIterator & operator=( const AlignataIterator & other ) {
      if (*this != other) {
	/* mIterator = other.mIterator; does not work, I guess, because of shallow copying */
	delete mIterator;
	mIterator = other.mIterator->getClone();
      }
      
      return *this;
    }

    /** comparison operator */
    inline bool operator==( const AlignataIterator & other) {   
      return ( (*mIterator).getPointer() == (*other.mIterator).getPointer() );
    }
    
    /** comparison operator */
    inline bool operator!=( const AlignataIterator & other) {
      return ( (*mIterator).getPointer() != (*other.mIterator).getPointer() );
    }
      
    /** indirection operator */
    inline ResiduePAIR * operator->() const { return &(mIterator->getReference()); }

    /** dereference operator */
    inline ResiduePAIR & operator*() const { return mIterator->getReference(); }
 
    /** prefix ++ */
    inline AlignataIterator & operator++() { mIterator->next(); return *this; }
 
    /** postfix ++ */
    inline AlignataIterator operator++(int) { AlignataIterator tmp = *this; mIterator->next(); return tmp; }
 
    /** prefix -- */
    inline AlignataIterator & operator--() { mIterator->previous(); return *this; }
 
    /** postfix -- */
    inline AlignataIterator  operator--(int) { 
      AlignataIterator tmp = *this; mIterator->previous(); return tmp; }
    
 private:
    Alignata::Iterator * mIterator;
};



}

#endif /* _ALIGNATA_H */

