/*
  alignlib - a library for aligning protein sequences

  $Id: AlignataIterator.cpp,v 1.4 2004/06/02 12:11:37 aheger Exp $

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
#include "alignlib.h"
#include "AlignlibDebug.h"
#include "AlignataIterator.h"

#ifdef WITH_DMALLOC
#include <dmalloc.h>
#endif

using namespace std;


namespace alignlib {

//------------------------------------< constructors and destructors >-----
AlignataConstIterator::AlignataConstIterator( Alignata::ConstIterator * it ) /* : mIterator(it) */ 
{
  debug_func_cerr(5);

  mIterator = it;
}

AlignataConstIterator::AlignataConstIterator( const AlignataConstIterator& src) : mIterator(0) 
{
  debug_func_cerr(5);


  mIterator = src.mIterator->getClone();
}

AlignataConstIterator::~AlignataConstIterator( ) {
    delete mIterator;
}

AlignataConstIterator & AlignataConstIterator::operator=( const AlignataConstIterator & other ) 
{
  debug_func_cerr(5);


  if (*this != other) {
    /* mIterator = other.mIterator; does not work, I guess, because of shallow copying */
    delete mIterator;
    mIterator = other.mIterator->getClone();
  }
  
  return *this;
}


const ResiduePAIR & AlignataConstIterator::operator*() const { return mIterator->getReference(); }


const ResiduePAIR * AlignataConstIterator::operator->() const { return &(mIterator->getReference()); }

  /* How to compare two iterators?
     1.Comparing the two mIterators would result in simply comparing the
     two memory locations of the iterators, which would not be valid. 
     2. use the operator == of ResiduePAIR. This seems to be the best choice, two iterators being equivalent, when they
     point to the same pair. Comparing the two residue pairs is straight-forward, however accessing them via reference breaks, 
     when the iterator points to nothing (for example, the iterator obtained by end()). Also, this is logically not correct.
     Two iterators working on different containers, but being on the same pair, would then be equal.
     3. Compare the memory location of the pairs. It is assumed, that end() will point to NULL
  */     
bool AlignataConstIterator::operator==( const AlignataConstIterator & other) { 
  return ( (*mIterator).getPointer() == (*other.mIterator).getPointer() ); 
}

bool AlignataConstIterator::operator!=( const AlignataConstIterator & other) { 
  return ( (*mIterator).getPointer() != (*other.mIterator).getPointer() ); 
}

AlignataConstIterator & AlignataConstIterator::operator++() { mIterator->next(); return *this; };
 
/** postfix ++ */
AlignataConstIterator AlignataConstIterator::operator++(int) { AlignataConstIterator tmp = *this; mIterator->next(); return tmp; }
 
/** prefix -- */
AlignataConstIterator & AlignataConstIterator::operator--() { mIterator->previous(); return *this; };
 
/** postfix -- */
AlignataConstIterator AlignataConstIterator::operator--(int) { AlignataConstIterator tmp = *this; mIterator->previous(); return tmp; }


//--------------------------------------> AlignataIterator <-----------------------------------------------------------------------


//------------------------------------< constructors and destructors >-----
AlignataIterator::AlignataIterator( Alignata::Iterator * it ) /* : mIterator(it) */ 
{
  debug_func_cerr(5);

  mIterator = it;
}

AlignataIterator::AlignataIterator( const AlignataIterator& src) : mIterator(0) 
{
  debug_func_cerr(5);


  mIterator = src.mIterator->getClone();
}

AlignataIterator::~AlignataIterator( ) {
    delete mIterator;
}

AlignataIterator & AlignataIterator::operator=( const AlignataIterator & other ) 
{
  debug_func_cerr(5);


  if (*this != other) {
    /* mIterator = other.mIterator; does not work, I guess, because of shallow copying */
    delete mIterator;
    mIterator = other.mIterator->getClone();
  }
  
  return *this;
}


ResiduePAIR & AlignataIterator::operator*() const { return mIterator->getReference(); }


ResiduePAIR * AlignataIterator::operator->() const { return &(mIterator->getReference()); }

  /* How to compare two iterators?
     1.Comparing the two mIterators would result in simply comparing the
     two memory locations of the iterators, which would not be valid. 
     2. use the operator == of ResiduePAIR. This seems to be the best choice, two iterators being equivalent, when they
     point to the same pair. Comparing the two residue pairs is straight-forward, however accessing them via reference breaks, 
     when the iterator points to nothing (for example, the iterator obtained by end()). Also, this is logically not correct.
     Two iterators working on different containers, but being on the same pair, would then be equal.
     3. Compare the memory location of the pairs. It is assumed, that end() will point to NULL
  */     
bool AlignataIterator::operator==( const AlignataIterator & other) { 
  return ( (*mIterator).getPointer() == (*other.mIterator).getPointer() ); 
}

bool AlignataIterator::operator!=( const AlignataIterator & other) { 
  return ( (*mIterator).getPointer() != (*other.mIterator).getPointer() ); 
}

AlignataIterator & AlignataIterator::operator++() { mIterator->next(); return *this; };
 
/** postfix ++ */
AlignataIterator AlignataIterator::operator++(int) { AlignataIterator tmp = *this; mIterator->next(); return tmp; }
 
/** prefix -- */
AlignataIterator & AlignataIterator::operator--() { mIterator->previous(); return *this; };
 
/** postfix -- */
AlignataIterator AlignataIterator::operator--(int) { AlignataIterator tmp = *this; mIterator->previous(); return tmp; }


} // namespace alignlib
