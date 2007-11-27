//--------------------------------------------------------------------------------
// Project LibAlign
//
// Copyright (C) 2000 Andreas Heger All rights reserved
//
// Author: Andreas Heger <heger@ebi.ac.uk>
//
// $Id: AlignedBlocks.cpp,v 1.1.1.1 2002/07/08 21:20:17 heger Exp $
//--------------------------------------------------------------------------------    


#include <iostream>
#include <iomanip>
#include <iterator>
#include <cstring>
#include <string>

#include "AlignlibDebug.h"
#include "Alignata.h"
#include "AlignataIterator.h"
#include "AlignedBlocks.h"
using namespace std;

namespace alignlib 
{

//---------------------------------------------------------< constructors and destructors >--------------------------------------
AlignedBlocks::AlignedBlocks (const Alignata * src) : 
	mRowFrom(NO_POS), mRowTo(NO_POS), mColFrom(NO_POS), mColTo(NO_POS)
	{
	if (src != NULL)
		fill( src );
	}


AlignedBlocks::~AlignedBlocks () 
{
	mRowStarts.clear();
	mColStarts.clear();
	mBlockSizes.clear();
}

AlignedBlocks::AlignedBlocks (const AlignedBlocks & src ) 
{
}

//---------------------------------------------------------------
void AlignedBlocks::fill( const Alignata * src)
{
	debug_func_cerr(5);

	mRowStarts.clear();
	mColStarts.clear();
	mBlockSizes.clear();

	mRowFrom = src->getRowFrom();
	mRowTo = src->getRowTo();
	mColFrom = src->getColFrom();
	mColTo = src->getColTo();

	// sanity checks
	if (src->isEmpty()) return;

	AlignataConstIterator it(src->begin());
	AlignataConstIterator it_end(src->end());

	Position last_col = it->mCol; 
	Position last_row = it->mRow; 

	Position d_row, d_col;

	// start iteration at col_from + 1
	mRowStarts.push_back( it->mRow - mRowFrom );
	mColStarts.push_back( it->mCol - mColFrom );
	Position block_size = 1;

	++it;

	for (; it != it_end; ++it)
	{
		Position current_row = it->mRow;
		Position current_col = it->mCol;

		if ( (current_row - last_row) > 1 || (current_col - last_col) > 1) 
		{
			mBlockSizes.push_back( block_size );
			mRowStarts.push_back( current_row - mRowFrom );
			mColStarts.push_back( current_col - mColFrom );
			block_size = 0;
		}       
		++ block_size;
		last_row = current_row;
		last_col = current_col;
	}
	mBlockSizes.push_back( block_size );
}

//--------------------------------------------------------------------------------------------------------------------------------
Alignata * AlignedBlocks::copy( Alignata * dest ) const 
{
	debug_func_cerr(5);

	for (int x = 0; x < mRowStarts.size(); ++x)
	{
		Position row = mRowStarts[x] + mRowFrom; 
		Position col = mColStarts[x] + mColFrom;		
		for (int l = 0; l < mBlockSizes[x]; ++l, ++row, ++col)
			dest->addPair( row, col );
	}
}

//--------------------------------------------------------------------------------------------------------------------------------
std::ostream & operator<< (std::ostream & output, const AlignedBlocks & src) 
{
	output << src.mRowFrom << "\t" << src.mRowTo << "\t" << src.mColFrom << "\t" << src.mColTo << "\t";
	std::copy( src.mRowStarts.begin(), src.mRowStarts.end(), std::ostream_iterator<Position>(output, ","));
	output << "\t";
	std::copy( src.mColStarts.begin(), src.mColStarts.end(), std::ostream_iterator<Position>(output, ","));
	output << "\t";
	std::copy( src.mBlockSizes.begin(), src.mBlockSizes.end(), std::ostream_iterator<Position>(output, ","));
	return output;
}

//--------------------------------------------------------------------------------------------------------------------------------
inline void parseList( std::istream & input, PositionVector & dest )
{
	std::string delimiters(",");
	std::string str;
	input >> str;

	string::size_type last_pos = str.find_first_not_of( delimiters, 0);
	string::size_type pos     = str.find_first_of(delimiters, last_pos);

	while (string::npos != pos || string::npos != last_pos)
	{
		dest.push_back(atoi(str.substr(last_pos, pos - last_pos).c_str()));
		last_pos = str.find_first_not_of(delimiters, pos);
		pos = str.find_first_of(delimiters, last_pos);
	}	
}

std::istream & operator>> (std::istream & input, AlignedBlocks & dest) 
{
	input >> dest.mRowFrom >> dest.mRowTo >> dest.mColFrom >> dest.mColTo;

	parseList( input, dest.mRowStarts );
	parseList( input, dest.mColStarts );
	parseList( input, dest.mBlockSizes );
	return input;
}


} // namespace alignlib
