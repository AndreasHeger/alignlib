/*
  alignlib - a library for aligning protein sequences

  $Id: AlignmentIterator.cpp,v 1.4 2004/06/02 12:11:37 aheger Exp $

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


#include <iostream> 
#include <algorithm>
#include <set>
#include "alignlib_fwd.h"
#include "alignlib_interfaces.h"
#include "AlignlibDebug.h"
#include "AlignmentIterator.h"

#ifdef WITH_DMALLOC
#include <dmalloc.h>
#endif

using namespace std;


namespace alignlib {

//------------------------------------< constructors and destructors >-----
AlignmentConstIterator::AlignmentConstIterator( Alignment::ConstIterator * it ) /* : mIterator(it) */ 
{
  debug_func_cerr(5);

  mIterator = it;
}

AlignmentConstIterator::AlignmentConstIterator( const AlignmentConstIterator& src) : mIterator(0) 
{
  debug_func_cerr(5);


  mIterator = src.mIterator->getClone();
}

AlignmentConstIterator::~AlignmentConstIterator( ) {
    delete mIterator;
}

AlignmentConstIterator & AlignmentConstIterator::operator=( const AlignmentConstIterator & other ) 
{
  debug_func_cerr(5);


  if (*this != other) {
    /* mIterator = other.mIterator; does not work, I guess, because of shallow copying */
    delete mIterator;
    mIterator = other.mIterator->getClone();
  }
  
  return *this;
}


const ResiduePair & AlignmentConstIterator::operator*() const { return mIterator->getReference(); }


const ResiduePair * AlignmentConstIterator::operator->() const { return &(mIterator->getReference()); }

  /* How to compare two iterators?
     1.Comparing the two mIterators would result in simply comparing the
     two memory locations of the iterators, which would not be valid. 
     2. use the operator == of ResiduePair. This seems to be the best choice, two iterators being equivalent, when they
     point to the same pair. Comparing the two residue pairs is straight-forward, however accessing them via reference breaks, 
     when the iterator points to nothing (for example, the iterator obtained by end()). Also, this is logically not correct.
     Two iterators working on different containers, but being on the same pair, would then be equal.
     3. Compare the memory location of the pairs. It is assumed, that end() will point to NULL
  */     
bool AlignmentConstIterator::operator==( const AlignmentConstIterator & other) { 
  return ( (*mIterator).getPointer() == (*other.mIterator).getPointer() ); 
}

bool AlignmentConstIterator::operator!=( const AlignmentConstIterator & other) { 
  return ( (*mIterator).getPointer() != (*other.mIterator).getPointer() ); 
}

AlignmentConstIterator & AlignmentConstIterator::operator++() { mIterator->next(); return *this; };
 
/** postfix ++ */
AlignmentConstIterator AlignmentConstIterator::operator++(int) { AlignmentConstIterator tmp = *this; mIterator->next(); return tmp; }
 
/** prefix -- */
AlignmentConstIterator & AlignmentConstIterator::operator--() { mIterator->previous(); return *this; };
 
/** postfix -- */
AlignmentConstIterator AlignmentConstIterator::operator--(int) { AlignmentConstIterator tmp = *this; mIterator->previous(); return tmp; }


//--------------------------------------> AlignmentIterator <-----------------------------------------------------------------------


//------------------------------------< constructors and destructors >-----
AlignmentIterator::AlignmentIterator( Alignment::Iterator * it ) /* : mIterator(it) */ 
{
  debug_func_cerr(5);

  mIterator = it;
}

AlignmentIterator::AlignmentIterator( const AlignmentIterator& src) : mIterator(0) 
{
  debug_func_cerr(5);


  mIterator = src.mIterator->getClone();
}

AlignmentIterator::~AlignmentIterator( ) {
    delete mIterator;
}

AlignmentIterator & AlignmentIterator::operator=( const AlignmentIterator & other ) 
{
  debug_func_cerr(5);


  if (*this != other) {
    /* mIterator = other.mIterator; does not work, I guess, because of shallow copying */
    delete mIterator;
    mIterator = other.mIterator->getClone();
  }
  
  return *this;
}


ResiduePair & AlignmentIterator::operator*() const { return mIterator->getReference(); }


ResiduePair * AlignmentIterator::operator->() const { return &(mIterator->getReference()); }

  /* How to compare two iterators?
     1.Comparing the two mIterators would result in simply comparing the
     two memory locations of the iterators, which would not be valid. 
     2. use the operator == of ResiduePair. This seems to be the best choice, two iterators being equivalent, when they
     point to the same pair. Comparing the two residue pairs is straight-forward, however accessing them via reference breaks, 
     when the iterator points to nothing (for example, the iterator obtained by end()). Also, this is logically not correct.
     Two iterators working on different containers, but being on the same pair, would then be equal.
     3. Compare the memory location of the pairs. It is assumed, that end() will point to NULL
  */     
bool AlignmentIterator::operator==( const AlignmentIterator & other) { 
  return ( (*mIterator).getPointer() == (*other.mIterator).getPointer() ); 
}

bool AlignmentIterator::operator!=( const AlignmentIterator & other) { 
  return ( (*mIterator).getPointer() != (*other.mIterator).getPointer() ); 
}

AlignmentIterator & AlignmentIterator::operator++() { mIterator->next(); return *this; };
 
/** postfix ++ */
AlignmentIterator AlignmentIterator::operator++(int) { AlignmentIterator tmp = *this; mIterator->next(); return tmp; }
 
/** prefix -- */
AlignmentIterator & AlignmentIterator::operator--() { mIterator->previous(); return *this; };
 
/** postfix -- */
AlignmentIterator AlignmentIterator::operator--(int) { AlignmentIterator tmp = *this; mIterator->previous(); return tmp; }


} // namespace alignlib
