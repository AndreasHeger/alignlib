/*
  alignlib - a library for aligning protein sequences

  $Id: ImplAlignmentSet.cpp,v 1.4 2004/06/02 12:11:37 aheger Exp $

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
#include "ImplAlignmentSorted.h"
#include "AlignmentIterator.h"

using namespace std;

namespace alignlib
{

//-- Comparator objects - these define the sort order of residues within an alignment

/** residues are sorted by row only */
struct ComparatorRow
{
  bool operator()( const ResiduePAIR * x, const ResiduePAIR * y) const 
  { 
    return x->mRow < y->mRow; 
  } 
};

/** residues are sorted by col only */
struct ComparatorCol
{
  bool operator()( const ResiduePAIR * x, const ResiduePAIR * y) const 
  { 
    return x->mCol < y->mCol; 
  } 
};

/** residues are sorted by row and then by column */
struct ComparatorRowCol
{
  bool operator() ( const ResiduePAIR * x, const ResiduePAIR * y) const { 
    if (x->mRow < y->mRow) return 1;
    if (x->mRow > y->mRow) return 0;
    if (x->mCol < y->mCol) 
      return 1;
    else
      return 0;
  }
};  

 /** residues are sorted by diagonal and then by column */
struct ComparatorDiagonalCol {
     bool operator() ( const ResiduePAIR * x, const ResiduePAIR * y) const { 
	  Diagonal d1 = x->mCol - x->mRow;
	  Diagonal d2 = y->mCol - y->mRow;
	  if (d1 < d2) return 1;
	  if (d1 > d2) return 0;
	  if (x->mCol < y->mCol) 
	      return 1;
	  else
	      return 0;
     }
 };  


typedef std::set<ResiduePAIR *, ComparatorRow> AlignmentSetRow;
typedef std::set<ResiduePAIR *, ComparatorCol> AlignmentSetCol;
typedef std::set<ResiduePAIR *, ComparatorRowCol> AlignmentSetRowCol;
typedef std::set<ResiduePAIR *, ComparatorDiagonalCol> AlignmentSetDiagonalCol;  

//------------------------------factory functions -----------------------------
Alignment * makeAlignmentSet()
{
	return new ImplAlignmentSorted< AlignmentSetRow >();
}

Alignment * makeAlignmentSetCol()
{
	return new ImplAlignmentSorted< AlignmentSetCol >();
}

Alignment * makeAlignmentHash()
{
	return new ImplAlignmentSorted< AlignmentSetRowCol >();
}

Alignment * makeAlignmentHashDiagonal()
{
	return new ImplAlignmentSorted< AlignmentSetDiagonalCol >();
}


//------------------------------------< constructors and destructors >-----
template <class T>
ImplAlignmentSorted<T>::ImplAlignmentSorted() : ImplAlignment() 
{
	debug_func_cerr(5);

}

template <class T>
ImplAlignmentSorted<T>::ImplAlignmentSorted( const ImplAlignmentSorted& src) : ImplAlignment( src ) 
{
	debug_func_cerr(5);


	clearContainer();

	PairIterator it(src.mPairs.begin()), it_end(src.mPairs.end());
	for (;it != it_end; ++it) 
		mPairs.insert( new ResiduePAIR(**it) );
}

template <class T>  
ImplAlignmentSorted<T>::~ImplAlignmentSorted( ) 
{
	debug_func_cerr(5);

	clear();
}

//------------------------------------------------------------------------------------------------------------
template <class T>
ImplAlignmentSorted<T> * ImplAlignmentSorted<T>::getNew() const
{
	return new ImplAlignmentSorted<T>();
}

template <class T>
ImplAlignmentSorted<T> * ImplAlignmentSorted<T>::getClone() const
{
	return new ImplAlignmentSorted( *this );
}

//------------------------------------------------------------------------------------------------------------
template <class T>  
void ImplAlignmentSorted<T>::clearContainer()
{ 
	PairIterator it(mPairs.begin()), it_end(mPairs.end());
	for (;it != it_end; ++it) delete *it;
	mPairs.clear(); 
}

//------------------------------------------------------------------------------------------------------------
template <class T>  
void ImplAlignmentSorted<T>::clear() { 
	ImplAlignment::clear();
	clearContainer();
}

//-----------------------------------------------------------------------------------------------------------   
template <class T>
AlignmentConstIterator ImplAlignmentSorted<T>::begin() const
{ 
	return AlignmentConstIterator( new ImplAlignmentSorted_ConstIterator(mPairs.begin()));
}

template <class T>  
AlignmentConstIterator ImplAlignmentSorted<T>::end()   const
{ 
	return AlignmentConstIterator( new ImplAlignmentSorted_ConstIterator(mPairs.end())); 
}

template <class T>
AlignmentIterator ImplAlignmentSorted<T>::begin()
{ 
	return AlignmentIterator( new ImplAlignmentSorted_Iterator(mPairs.begin()));
}

template <class T>
AlignmentIterator ImplAlignmentSorted<T>::end()
{ 
	return AlignmentIterator( new ImplAlignmentSorted_Iterator(mPairs.end())); 
}


//----------------> accessors <------------------------------------------------------------------------------

template <class T>
ResiduePAIR ImplAlignmentSorted<T>::front() const { return **(mPairs.begin()); }
template <class T>  
ResiduePAIR ImplAlignmentSorted<T>::back()  const { return **(mPairs.rbegin()); }

template <class T>
void ImplAlignmentSorted<T>::addPair( ResiduePAIR * new_pair ) 
{
	debug_func_cerr(5);

	ImplAlignment::addPair( new_pair );
	if (mPairs.find( new_pair) != mPairs.end()) 
	{
		delete new_pair;
	} else 
	{
		setChangedLength(); 
		mPairs.insert( new_pair ); 
	} 
} 

//----------------------------------------------------------------------------------------------------------
/** retrieves a pair of residues from the alignment */
template <class T>  
ResiduePAIR ImplAlignmentSorted<T>::getPair( const ResiduePAIR & p) const 
{
	PairIterator it = mPairs.find( const_cast< ResiduePAIR* > (&p) );
	if (it != mPairs.end())
		return **it;
	else
		return ResiduePAIR();
} 

//----------------------------------------------------------------------------------------------------------
/** remove a pair from an alignment */
template <class T>  
void ImplAlignmentSorted<T>::removePair( const ResiduePAIR & p ) 
{
	debug_func_cerr(5);

	PairIterator it(mPairs.find( const_cast< ResiduePAIR*> (&p)) );

	if (it != mPairs.end()) 
	{
		setChangedLength(); 
		delete *it;
		mPairs.erase(it);
	}
} 

//----------------------------------------------------------------------------------------
/** This is non-generic routine. Since pairs are accessed by row, this is quite quick.
 */
template <class T>  
void ImplAlignmentSorted<T>::removeRowRegion( Position from, Position to) 
{

	for (Position pos = from; pos < to; pos++) 
	{
		ResiduePAIR p(pos, NO_POS, 0);
		PairIterator it(mPairs.find( &p ));

		if (it != mPairs.end()) {
			setChangedLength(); 
			delete *it;
			mPairs.erase(it);
		}
	}

}

//----------------------------------------------------------------------------------------
/** This is a generic routine. It creates a new alignment by making a copy of the old one */
template <class T>  
void ImplAlignmentSorted<T>::removeColRegion( Position from, Position to) 
{

	PairIterator it(mPairs.begin()), it_end(mPairs.end());

	// Valgrind did not like the iterator being deleted,
	// thus this complicated loop structure. Did not complain 
	while (it != it_end) 
	{
		if ( (*it)->mCol >= from && (*it)->mCol < to) 
		{
			delete *it;
			PairIterator it2 = it;
			++it;              
			mPairs.erase(it2);
		}
		else  
			++it;           
	}

	setChangedLength();

}

} // namespace alignlib
