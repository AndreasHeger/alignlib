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
#include "AlignmentFormat.h"

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
	mIsAligned.resize( mLength, false);
	AlignmentIterator it( map_mali2sequence->begin() ), end( map_mali2sequence->end());
	for( ; it != end; ++it) mIsAligned[it->mRow] = true;
	return;
}

//-----------------------------------------------------------------------------------------------------------
void ImplMultAlignment::buildAligned()
{
	debug_func_cerr(5);
	mIsAligned.clear();
	mIsAligned.resize( mLength, false);
	for (int x = 0; x < mRows.size(); ++ x)
	{
		AlignmentIterator it( mRows[x]->begin() ), end( mRows[x]->end());
		for( ; it != end; ++it) mIsAligned[it->mRow] = true;
	}
	return;
}

//------------------------------------------------------------------------------------
/** Add a full multiple alignment to the another alignment.
 */
void ImplMultAlignment::add(
		const HMultAlignment & other,
		const HAlignment & map_this2other )
{
	debug_func_cerr(5);

	// do not add empty mali
	if (other->isEmpty()) return;

	for (int x = 0; x < other->getNumSequences(); ++x)
	{
		HAlignment new_map_mali2sequence(other->getRow(x)->getNew());

		combineAlignment(
				new_map_mali2sequence,
				map_this2other,
				other->getRow(x),
				CR);

		mRows.push_back( new_map_mali2sequence );
	}

	mLength = std::max( mLength, map_this2other->getRowTo() );
	buildAligned();
}

//------------------------------------------------------------------------------------
/** Add a full multiple alignment to the another alignment.
 */
void ImplMultAlignment::add(
		const HMultAlignment & other,
		const HAlignment & map_this2new,
		const HAlignment & map_other2new )
{
	debug_func_cerr(5);

	// do not add empty mali
	if (other->isEmpty()) return;

	// map this alignment
	for (int x = 0; x < getNumSequences(); ++x)
	{
		mRows[x]->map( map_this2new, RR );
	}

	for (int x = 0; x < other->getNumSequences(); ++x)
	{
		HAlignment new_map_mali2sequence(other->getRow(x)->getClone());
		new_map_mali2sequence->map( map_other2new, RR);
		mRows.push_back( new_map_mali2sequence );
	}

	mLength = std::max( map_this2new->getColTo(), map_other2new->getColTo() );
	buildAligned();
}

//------------------------------------------------------------------------------------
/* add single entry to *this multiple alignment given an alignment.
 */
void ImplMultAlignment::add(
		const HAlignment & map_mali2sequence )
{
	debug_func_cerr(5);
	mRows.push_back( map_mali2sequence->getClone() );
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


//---------------------------------------------------------------------------------------
void ImplMultAlignment::expand( const HAlignandumVector & sequences )
{

	debug_func_cerr(5);

	if (isEmpty()) return;

	bool insert_termini = false;
	if (sequences->size() != 0)
	{
		if (sequences->size() != getNumSequences())
			throw AlignlibException( "ImplMultAlignment.cpp: number of sequences given does not match number of sequences in MultAlignment");
		insert_termini = true;
	}

	Position mali_length = 0;

	// find number of aligned columns in mali
	for (unsigned int x = 0; x < mRows.size(); ++x)
		mali_length = std::max( mali_length, mRows[x]->getRowTo());

	// find total/maximum insertions before a given mali column
	// row: position in mali
	// col: position in sequence
	std::vector<int> gaps(mali_length + 1, 0);

	for (unsigned int x = 0; x < mRows.size(); ++x)
	{
		HAlignment map_mali2row = mRows[x];

		Position last_col = map_mali2row->getColFrom();

		if (insert_termini)
			gaps[0] += last_col;

		for (Position row = map_mali2row->getRowFrom() + 1; row < map_mali2row->getRowTo(); ++row)
		{
			Position col = map_mali2row->mapRowToCol(row);
			if (col != NO_POS)
			{
				gaps[row] += col - last_col - 1;
			}
			last_col = col;
		}

		if (insert_termini)
			gaps[mali_length] += (*sequences)[x]->getLength() - map_mali2row->getColTo();

	}

	debug_cerr( 5, "length=" << mali_length << " insert_termini=" << insert_termini);

#ifdef DEBUG
	for (unsigned int x = 0; x < gaps.size(); ++x)
		debug_cerr( 5, "col=" << x << " gaps=" << gaps[x]);
#endif

	// build map of aligned columns to output columns in mali
	// record gaps before each position
	HAlignment map_mali_old2new = makeAlignmentVector();
	{
		Position y = 0;
		for (Position x = 0; x < mali_length; ++x)
		{
			y += gaps[x];
			map_mali_old2new->addPair( x, y++, 0 );
		}
	}

	debug_cerr( 5, "map_mali_old2new\n" << *map_mali_old2new );

	// remap each row to the new mali
	mLength = 0;

	std::vector<int>used_gaps(mali_length + 1, 0);

	for (unsigned int x = 0; x < mRows.size(); ++x)
	{
		HAlignment old_map_mali2row = mRows[x];
		HAlignment new_map_mali2row = old_map_mali2row->getNew();

		// build new alignment by mapping the existing aligned columns
		combineAlignment( new_map_mali2row,
						  map_mali_old2new,
						  old_map_mali2row,
						  RR);

		debug_cerr( 5, "map_mali2row after mapping aligned columns=\n" << *new_map_mali2row );

		// add residues for unaligned positions
		// insert before start
		if (insert_termini)
		{
			unsigned int u = used_gaps[0];
			Position col = old_map_mali2row->getColFrom();
			Position s = 0;
			while (s < col)
			{
				assert( new_map_mali2row->mapRowToCol(u) == NO_POS);
				new_map_mali2row->addPair( u++, s++, 0);
			}
			used_gaps[0] = u;
		}

		// insert gaps between aligned positions
		Position last_col = old_map_mali2row->getColFrom();
		for (Position row = old_map_mali2row->getRowFrom() + 1; row < old_map_mali2row->getRowTo(); ++row)
		{
			Position col = old_map_mali2row->mapRowToCol(row);

			if (col != NO_POS)
			{
				unsigned int u = map_mali_old2new->mapRowToCol(row) - gaps[row] + used_gaps[row];
				unsigned int d = col - last_col - 1;
				while (col - last_col - 1 > 0)
				{
					new_map_mali2row->addPair( u++, ++last_col, 0);
				}
				used_gaps[row] += d;
				last_col = col;
			}
		}

		if (insert_termini)
		{
			Position end = map_mali_old2new->getColTo();
			Position residues = (*sequences)[x]->getLength() - old_map_mali2row->getColTo();
			Position start = end + used_gaps[ mali_length ];
			debug_cerr( 5, "adding terminal residues:"
					<< " end=" << end
					<< " start=" << start
					<< " to=" << start + residues
					<< " diag=" << old_map_mali2row->getColTo() - start);
			new_map_mali2row->addDiagonal(
					start,
					start + residues,
					old_map_mali2row->getColTo() - start );
			used_gaps[mali_length] += residues;
		}

		debug_cerr( 5, "map_mali2row after mapping unaligned columns=\n" << *new_map_mali2row );

		mRows[x] = new_map_mali2row;
		mLength = std::max( mLength, new_map_mali2row->getRowTo() );
	}


	// by definition all columns will be aligned
	mIsAligned.clear();
	mIsAligned.resize( mLength, true);
}


//---------------------------------------------------------< Input/Output routines >--------

HPositionMatrix ImplMultAlignment::getPositionMatrix( const bool & transpose ) const
{
	debug_func_cerr(5);
	PositionMatrix * matrix;
	if (transpose)
	{
		matrix = new PositionMatrix( getLength(), getNumSequences(), NO_POS);
		for( int x = 0; x < mRows.size(); ++x)
		{
			AlignmentIterator it(mRows[x]->begin()), end(mRows[x]->end());
			for( ; it!= end; ++it)
				matrix->setValue( it->mRow, x, it->mCol );
		}
	}
	else
	{
		matrix = new PositionMatrix( getNumSequences(), getLength(), NO_POS);
		for( int x = 0; x < mRows.size(); ++x)
		{
			AlignmentIterator it(mRows[x]->begin()), end(mRows[x]->end());
			for( ; it!= end; ++it)
				matrix->setValue( x, it->mRow, it->mCol );
		}
	}
	return HPositionMatrix( matrix );
}

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


