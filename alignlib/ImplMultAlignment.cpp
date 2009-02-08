/*
  alignlib - a library for aligning protein sequences

  $Id: ImplMultAlignment.cpp,v 1.6 2004/03/19 18:23:41 aheger Exp $

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
#include <fstream>
#include <iomanip>
#include <vector>
#include <cassert>
#include "alignlib_fwd.h"
#include "alignlib_interfaces.h"
#include "AlignlibDebug.h"

#include "HelpersAlignandum.h"
#include "ImplMultAlignment.h"
#include "HelpersMultAlignment.h"
#include "HelpersEncoder.h"
#include "AlignlibException.h"
#include "Alignatum.h"
#include "Alignandum.h"
#include "Alignment.h"
#include "HelpersAlignment.h"
#include "AlignmentIterator.h"
#include "AlignlibException.h"
#include "HelpersAlignment.h"

using namespace std;

namespace alignlib
{

/* factory functions */

/** create an empty multiple alignment */
HMultAlignment makeMultAlignment()
{
	return HMultAlignment( new ImplMultAlignment() );
}


//---------------------------------------------------------< constructors and destructors >--------------------------------------
ImplMultAlignment::ImplMultAlignment () :
	mLength(0)
	{
	}

//--------------------------------------------------------------------------------------------------------------
ImplMultAlignment::~ImplMultAlignment ()
{
	freeMemory();
}

//--------------------------------------------------------------------------------------------------------------
ImplMultAlignment::ImplMultAlignment (const ImplMultAlignment & src )
:
	mLength (src.mLength)
{

	// clear old entries
	freeMemory();

	// add clones of the new entries
	for (unsigned int row = 0; row < src.mRows.size(); row++)
		add( src.mRows[row]->getClone() );

	mIsAligned.clear();
	std::copy(src.mIsAligned.begin(), src.mIsAligned.end(), std::back_inserter( mIsAligned ));
}

//--------------------------------------------------------------------------------------------------------------
void ImplMultAlignment::freeMemory()
{
	debug_func_cerr(5);
	mRows.clear();
	mIsAligned.clear();
}

//--------------------------------------------------------------------------------------------------------------
Position ImplMultAlignment::getLength() const
{
	return mLength;
}

//--------------------------------------------------------------------------------------------------------------
int ImplMultAlignment::getNumSequences() const
{
	return mRows.size();
}

//-----------------------------------------------------------------------------------------------------------
const HAlignment ImplMultAlignment::operator[]( int row ) const
{
	if (isEmpty())
		throw AlignlibException("In ImplMultAlignment.cpp: alignment is empty");

	if (row < 0 || row >= mRows.size() )
		throw AlignlibException("In ImplMultAlignment.cpp: out-of-range access");

	return mRows[row];
}

//-----------------------------------------------------------------------------------------------------------
const HAlignment ImplMultAlignment::getRow( int row ) const
{
	debug_func_cerr(5);
	if (isEmpty())
		throw AlignlibException("In ImplMultAlignment.cpp: alignment is empty");
	if (row < 0 || row >= mRows.size() )
		throw AlignlibException("In ImplMultAlignment.cpp: out-of-range access");

	return mRows[row];
}

//-----------------------------------------------------------------------------------------------------------
void ImplMultAlignment::clear()
{
	debug_func_cerr(5);
	freeMemory();
	mLength = 0;
}

//-----------------------------------------------------------------------------------------------------------
void ImplMultAlignment::eraseRow( int row )
{
	debug_func_cerr(5);
	if (isEmpty())
		throw AlignlibException("In ImplMultAlignment.cpp: alignment is empty");
	if (row < 0 || row >= mRows.size() )
		throw AlignlibException("In ImplMultAlignment.cpp: out-of-range access");

	mRows.erase( mRows.begin() + row );
	if (mRows.size() == 0)
		mLength = 0;
}

//-----------------------------------------------------------------------------------------------------------
bool ImplMultAlignment::isAligned( const Position & col )
{
	debug_func_cerr(5);
	if (col < 0 || col >= getLength())
		throw AlignlibException("In ImplMultAlignment.cpp: out-of-range access");
	return mIsAligned[ col ];
}


//-----------------------------------------------------------------------------------------------------------
void ImplMultAlignment::updateAligned(
		const HAlignment & map_mali2sequence )
{
	debug_func_cerr(5);
	mIsAligned.clear();
	mIsAligned.resize( mLength, false);
	AlignmentIterator it( map_mali2sequence->begin() ), end( map_mali2sequence->end());
	for( ; it != end; ++it) mIsAligned[it->mRow] = true;
	return;
}

//------------------------------------------------------------------------------------
/** Add a full multiple alignment to the another alignment.
 */
void ImplMultAlignment::add(
		const HMultAlignment & src,
		const HAlignment & map_mali2sequence )
{
	debug_func_cerr(5);

	// do not add empty mali
	if (src->isEmpty()) return;

	for (int x = 0; x < src->getNumSequences(); ++x)
	{
		HAlignment new_map_mali2sequence(src->getRow(x)->getNew());

		combineAlignment(
				new_map_mali2sequence,
				map_mali2sequence,
				src->getRow(x),
				CR);

		mRows.push_back( new_map_mali2sequence );
	}

	mLength = std::max( mLength, map_mali2sequence->getRowTo() );
	updateAligned( map_mali2sequence );
}

//------------------------------------------------------------------------------------
/* add single entry to *this multiple alignment given an alignment.
 */
void ImplMultAlignment::add(
		const HAlignment & map_mali2sequence )
{
	debug_func_cerr(5);
	mRows.push_back( map_mali2sequence );
	mLength = std::max( mLength, map_mali2sequence->getRowTo() );
	updateAligned( map_mali2sequence );
}

//---------------------------------------------------------------------------------------
HMultAlignment ImplMultAlignment::getClone() const
{
	return HMultAlignment( new ImplMultAlignment( *this ) );
}

//---------------------------------------------------------------------------------------
HMultAlignment ImplMultAlignment::getNew() const
{
	return HMultAlignment( new ImplMultAlignment( ) );
}

//---------------------------------------------------------------------------------------
bool ImplMultAlignment::isEmpty() const
{
	return mRows.empty();
}

//---------------------------------------------------------< Input/Output routines >--------

//------------------------------------------------------------------------------------
void ImplMultAlignment::write( std::ostream & output ) const
{
	debug_func_cerr(5);

	for (unsigned int row = 0; row < mRows.size(); ++row)
	{
		mRows[row]->write( output );
		output << std::endl;
	}
	Position l = getLength();
	for (unsigned int col = 0; col < l; ++col)
		output << mIsAligned[col];
}

} // namespace alignlib


