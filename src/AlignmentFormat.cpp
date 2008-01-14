//--------------------------------------------------------------------------------
// Project LibAlign
//
// Copyright (C) 2000 Andreas Heger All rights reserved
//
// Author: Andreas Heger <heger@ebi.ac.uk>
//
// $Id: AlignmentFormatBlocks.cpp,v 1.1.1.1 2002/07/08 21:20:17 heger Exp $
//--------------------------------------------------------------------------------    


#include <iostream>
#include <iomanip>
#include <iterator>
#include <cstring>
#include <string>

#include "alignlib.h"
#include "AlignlibDebug.h"
#include "AlignException.h"
#include "AlignmentFormatBlocks.h"

using namespace std;

namespace alignlib 
{

//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
AlignmentFormat::AlignmentFormat() : 
	mRowFrom(NO_POS), mRowTo(NO_POS), mColFrom(NO_POS), mColTo(NO_POS)
	{
	}

AlignmentFormat::AlignmentFormat( const AlignmentFormat & src) : 
	mRowFrom(src.mRowFrom), mRowTo(src.mRowTo), 
	mColFrom(src.mColFrom), mColTo(src.mColTo)
	{
	}

AlignmentFormat::~AlignmentFormat()
{
}

void AlignmentFormat::fill( const HAlignment & src )
{
	mRowFrom = src->getRowFrom();
	mRowTo = src->getRowTo();
	mColFrom = src->getColFrom();
	mColTo = src->getColTo();
}

void AlignmentFormat::copy( HAlignment & dest )
{
	dest->clear();
}

//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
AlignmentFormatBlocks::AlignmentFormatBlocks() : 
	AlignmentFormat()
	{
	}

AlignmentFormatBlocks::AlignmentFormatBlocks( const HAlignment & src) : 
	AlignmentFormat()
	{
	fill( src );
	}


AlignmentFormatBlocks::~AlignmentFormatBlocks () 
{
	mRowStarts.clear();
	mColStarts.clear();
	mBlockSizes.clear();
}

AlignmentFormatBlocks::AlignmentFormatBlocks (const AlignmentFormatBlocks & src ) :
	AlignmentFormat( src )
{
}

void AlignmentFormatBlocks::fill( const HAlignment & src)
{
	debug_func_cerr(5);

	AlignmentFormat::fill( src );
	
	mRowStarts.clear();
	mColStarts.clear();
	mBlockSizes.clear();

	// sanity checks
	if (src->isEmpty()) return;

	AlignmentIterator it(src->begin());
	AlignmentIterator it_end(src->end());

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
void AlignmentFormatBlocks::copy( HAlignment & dest ) const 
{
	debug_func_cerr(5);

	AlignmentFormat::copy( dest );
	
	for (int x = 0; x < mRowStarts.size(); ++x)
	{
		Position row = mRowStarts[x] + mRowFrom; 
		Position col = mColStarts[x] + mColFrom;		
		for (int l = 0; l < mBlockSizes[x]; ++l, ++row, ++col)
			dest->addPair( row, col );
	}
}

//--------------------------------------------------------------------------------------------------------------------------------
std::ostream & operator<< (std::ostream & output, const AlignmentFormatBlocks & src) 
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

std::istream & operator>> (std::istream & input, AlignmentFormatBlocks & dest) 
{
	input >> dest.mRowFrom >> dest.mRowTo >> dest.mColFrom >> dest.mColTo;

	parseList( input, dest.mRowStarts );
	parseList( input, dest.mColStarts );
	parseList( input, dest.mBlockSizes );
	return input;
}

//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
AlignmentFormatEmissions::AlignmentFormatEmissions() : 
	AlignmentFormat(),
	mRowAlignment(""),
	mColAlignment("")	
	{
	}

AlignmentFormatEmissions::AlignmentFormatEmissions( const HAlignment & src) : 
	AlignmentFormat()
	{
	fill( src );
	}


AlignmentFormatEmissions::~AlignmentFormatEmissions () 
{
}

AlignmentFormatEmissions::AlignmentFormatEmissions (const AlignmentFormatEmissions & src ) :
	AlignmentFormat( src ), 
	mRowAlignment( src.mRowAlignment ), 
	mColAlignment( src.mColAlignment )
{
}

void AlignmentFormatEmissions::fill( const HAlignment & src)
{
	debug_func_cerr(5);

	AlignmentFormat::fill( src );
	
	// sanity checks
	if (src->isEmpty()) return;

	AlignmentIterator it(src->begin());
	AlignmentIterator it_end(src->end());

	std::ostringstream os_row;
	std::ostringstream os_col;

	Position last_col = it->mCol; 
	Position len_col = 0;
	
	while (last_col < col_from) 
	{
		++it;
		last_col = it->mCol;
	}
	
	Position last_row = it->mRow; Position len_row = 0;

	Position d_row, d_col;
	Position current_row = 0, current_col = 0;

	// start iteration at col_from + 1
	++it; len_row++; len_col++;
	for (; it != it_end; ++it)
	{
		current_row = it->mRow;
		current_col = it->mCol;

		debug_cerr( 5, "current_row: " << current_row << " " << "current_col: " << current_col);

		if (current_col > col_to)
			break;

		if ((d_row = current_row - last_row - 1) > 0) 
		{
			os_col << "+" << len_col << "-" << d_row;

			len_col = 0;
			len_row += d_row;
		}

		if ((d_col = current_col - last_col - 1) > 0) 
		{
			os_row << "+" << len_row << "-" << d_col;
			len_row = 0;
			len_col += d_col;
		}

		last_row = current_row;
		last_col = current_col;
		len_row++; 
		len_col++;
	}

	// no ends necessary as with strstream,
	// ostringstream does it automatically.
	os_col << "+" << len_col;
	os_row << "+" << len_row;
	
	mRowAlignment = os_row;
	mColAlignment = os_col;
}

//--------------------------------------------------------------------------------------------------------------------------------
void AlignmentFormatEmissions::copy( HAlignment & dest ) const 
{
	debug_func_cerr(5);

	AlignmentFormat::copy( dest );
	
	if (mRowFrom == NO_POS || mColFrom == NO_POS)
		throw AlignException( "AlignmentFormat.cpp: alignment ranges not defined." );
	
	std::istringstream is_row( mRowAlignment.c_str() );   
	std::istringstream is_col( mColAlignment.c_str() );  

	Position row = mRowFrom;   
	Position col = mColFrom;   
	Position d_row = 0;   
	Position d_col = 0;  

	if (!(is_row >> d_row)) return;
	if (!(is_col >> d_col)) return;

	while (true) 
	{
		// cout << "0:" << d_row << " " << d_col << " " << row << " " << col << endl;
		// entry: d_row and d_col > 0
		// end: d_col or d_row or both are 0
		// emit pairs for aligned regions, i.e. d_row and d_col are positive
		while (d_row > 0 && d_col > 0) 
		{         
			dest->addPair( new ResiduePAIR (row, col) );         
			d_row--; 
			d_col--;         
			row++; 
			col++;       
		}  

		// cout << "1:" << d_row << " " << d_col << " " << row << " " << col << endl;

		// entry: d_row and d_col < 0
		// end: d_col or d_row or both are 0
		// simply skip. This region is a gap in both sequences
		if (d_row < 0 && d_col < 0) 
		{         
			if (d_row < d_col) 
			{
				d_row -= d_col;
				d_col = 0;
			} 
			else 
			{
				d_col -= d_row;
				d_row = 0;
			}
		}

		// cout << "2:" << d_row << " " << d_col << " " << row << " " << col << endl;

		// entry: d_row > 0, d_col < 0:
		// exit: d_row or d_col or both are 0
		// emitting characters only from row
		if (d_row > 0 && d_col < 0) 
		{
			if (d_row > -d_col ) 
			{
				d_row += d_col;
				row   -= d_col;		// emit
				d_col  = 0;
			} 
			else 
			{
				d_col += d_row;
				row   += d_row;		// emit
				d_row  = 0;
			}
		}

		// cout << "3:" << d_row << " " << d_col << " " << row << " " << col << endl;

		// entry: d_row < 0, d_col > 0:
		// exit: d_row or d_col or both are 0
		// emitting characters only from col
		if (d_col > 0 && d_row < 0) 
		{
			if (d_col > -d_row ) 
			{
				d_col += d_row;
				col   -= d_row;		// emit 
				d_row  = 0;
			} 
			else 
			{
				d_row += d_col;
				col   += d_col;		// emit
				d_col  = 0;
			}
		}

		// cout << "4:" << d_row << " " << d_col << " " << row << " " << col << endl;

		if (d_row == 0) 
			if (!(is_row >> d_row)) break;	// read new d_row

		if (d_col == 0) 
			if (!(is_col >> d_col)) break;	// read new d_col

	}	  

	return;	
}

//--------------------------------------------------------------------------------------------------------------------------------
std::ostream & operator<< (std::ostream & output, const AlignmentFormatEmissions & src) 
{
	output 
		<< src.mRowFrom << "\t" << src.mRowTo << "\t" << mRowAlignment << "\t" 
		<< src.mColFrom << "\t" << src.mColTo << "\t" << mColAlignment << "\t"
	return output;
}

std::istream & operator>> (std::istream & input, AlignmentFormatEmissions & dest) 
{
	input >> dest.mRowFrom >> dest.mRowTo >> dest.mRowAlignment >> dest.mColFrom >> dest.mColTo >> dest.mColAlignment;
	return input;
}

//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
AlignmentFormatExplicit::AlignmentFormatExplicit() : 
	AlignmentFormat(), 
	mRowAlignment(""),
	mColAlignment("")
	{
	}

AlignmentFormatExplicit::AlignmentFormatExplicit( 
		const HAlignment & src,
		const HAlignandum & row,
		const HAlignandum & col) : 
	AlignmentFormat()
	{
			fill( src, row, col );
	}


AlignmentFormatExplicit::~AlignmentFormatExplicit () 
{
}

AlignmentFormatExplicit::AlignmentFormatExplicit (const AlignmentFormatExplicit & src ) :
	AlignmentFormat( src ), 
	mRowAlignment( src.mRowAlignment ), 
	mColAlignment( src.mColAlignment )
{
}

void AlignmentFormatExplicit::fill( 
		const HAlignment & src,
		const HAlignandum & row,
		const HAlignandum & col )
{
	debug_func_cerr(5);

	AlignmentFormat::fill( src );
	
	// sanity checks
	if (src->isEmpty()) return;

	if (src->getRowTo() > row->getFullLength() )
		throw AlignException("alignment for row is out of bounds.");

	if (src->getColTo() > col->getFullLength() )	
		throw AlignException("alignment for col is out of bounds.");	

	HAlignment map_row2new = makeAlignmentVector();
	HAlignment map_col2new = makeAlignmentVector();

	fillAlignmentSummation( map_row2new, 
			map_col2new, 
			src, 
			true, true);

	HAlignatum row_alignatum = makeAlignatum( row, map_row2new);
	HAlignatum col_alignatum = makeAlignatum( col, map_col2new);

	mRowAlignment = row_alignatum->getString();
	mColAlignment = col_alignatum->getString();
	
	return;
}
	 	
//--------------------------------------------------------------------------------------------------------------------------------
void AlignmentFormatExplicit::copy( HAlignment & dest ) const 
{
	debug_func_cerr(5);

	AlignmentFormat::copy( dest );
	
	if (mRowFrom == NO_POS || mColFrom == NO_POS)
		throw AlignException( "AlignmentFormat.cpp: alignment ranges not defined." );

	char gap_char = getDefaultTranslator()->getGapChar();

	Position row = mRowFrom;   
	Position col = mColFrom;   

	for (unsigned int i = 0; i < row_ali.size(); i++) 
	{

		if (row_ali[i] != gap_char && col_ali[i] != gap_char) 
			dest->addPair( new ResiduePAIR (row, col) );         

		if (row_ali[i] != gap_char) 
			row++;

		if (col_ali[i] != gap_char) 
			col++;
	}

	return;	
}

//--------------------------------------------------------------------------------------------------------------------------------
std::ostream & operator<< (std::ostream & output, const AlignmentFormatExplicit & src) 
{
	output 
		<< src.mRowFrom << "\t" << src.mRowTo << "\t" << mRowAlignment << "\t" 
		<< src.mColFrom << "\t" << src.mColTo << "\t" << mColAlignment << "\t"
	return output;
}

std::istream & operator>> (std::istream & input, AlignmentFormatExplicit & dest) 
{
	input >> dest.mRowFrom >> dest.mRowTo >> dest.mRowAlignment >> dest.mColFrom >> dest.mColTo >> dest.mColAlignment;
	return input;
}

//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
AlignmentFormatDiagonals::AlignmentFormatDiagonals() : 
	AlignmentFormat(), mAlignment("")
	{
	}

AlignmentFormatDiagonals::AlignmentFormatDiagonals( 
		const HAlignment & src) : 
	AlignmentFormat()
	{
			fill( src);
	}


AlignmentFormatDiagonals::~AlignmentFormatDiagonals () 
{
}

AlignmentFormatDiagonals::AlignmentFormatDiagonals (const AlignmentFormatDiagonals & src ) :
	AlignmentFormat( src ), 
	mAlignment( src.mAlignment ), 
{
}

void AlignmentFormatDiagonals::fill( 
		const HAlignment & src,
		const bool reverse,
		const Position row_from,
		const Position row_to,
		const Position col_from,
		const Position col_to,
		const Diagonal diagonal_from,
		const Diagonal diagonal_to )
{
	debug_func_cerr(5);

	AlignmentFormat::fill( src );

	// sanity checks
	if (src->isEmpty()) return;

	// check parameters for filters and set them to sensible values
	if (col_from < src->getColFrom() || col_from == NO_POS)
		col_from = src->getColFrom();
	if (col_to >= src->getColTo() || col_to == NO_POS)
		col_to = src->getColTo();
	if (row_from < src->getRowFrom() || row_from == NO_POS)
		row_from = src->getRowFrom();
	if (row_to >= src->getRowTo() || row_to == NO_POS)
		row_to = src->getRowTo();

	if (diagonal_from > diagonal_to) 
	{
		diagonal_from = -MAX_DIAGONAL;
		diagonal_to   =  MAX_DIAGONAL;
	}

	// declare variables you need for iteration of the pairs
	AlignmentIterator it(src->begin());
	AlignmentIterator it_end(src->end());

	Diagonal last_diagonal = calculateDiagonal( *it );
	Diagonal this_diagonal = 0;

	Position this_row	 = it->mRow;
	Position this_col	 = it->mCol;
	Position last_row	 = this_row -1;
	Position emissions = 0;
	Position initial_gaps = 0;

	bool first = true;

	// now iterate over all pairs in the alignment
	for (;it!= it_end; ++it) {

		this_diagonal = calculateDiagonal(*it);
		this_row      = it->mRow;
		this_col      = it->mCol;

		debug_cerr( 5, "Pair:" << *it << std::endl            
				<< "last_diagonal=" << last_diagonal 
				<< " this_diagonal=" << this_diagonal
				<< " emissions=" << emissions 
				<< " this_row=" << this_row 
				<< " last_row=" << last_row );

		// apply filters
		if (this_col < col_from || this_col >= col_to) 
			continue;
		if (this_row < row_from || this_row >= row_to) 
			continue;
		if (this_diagonal < diagonal_from || this_diagonal > diagonal_to)
			continue;

		if (last_diagonal != this_diagonal || last_row >= this_row || first) 
		{

			if (!first) 
				output  << "+" << emissions << ";";

				// write last emission and switch to new diagonal
			if (this_diagonal < 0) 
				initial_gaps = this_col; 
			else
				initial_gaps = this_row;

			if (reverse)
				output << -this_diagonal << ":-" << initial_gaps;
			else
				output << this_diagonal << ":-" << initial_gaps;

			first = false;

			last_diagonal = this_diagonal;
			last_row = this_row;
			emissions = 1;
			continue;
		}

		// insert a gap
		if (last_row < this_row - 1) 
		{
			output << "+" << emissions << "-" << (this_row - last_row - 1);
			emissions = 0;
		}

		last_row = this_row;
		++emissions;

	}

	output << "+" << emissions;
	
	mAlignment = output;
	
	return;
}

//--------------------------------------------------------------------------------------------------------------------------------
void AlignmentFormatDiagonals::copy( HAlignment & dest ) const 
{
	debug_func_cerr(5);

	AlignmentFormat::copy( dest );
	
	if (mRowFrom == NO_POS || mColFrom == NO_POS)
		throw AlignException( "AlignmentFormat.cpp: alignment ranges not defined." );

	std::istringstream is_ali( ali.c_str() );   

	// set these/use these parameters to shift alignment
	Position row_from = mRowFrom;
	Position col_from = mColFrom;

	// row and col are the index of the next dot to be written.
	Position row = row_from;   
	Position col = col_from;   

	// read diagonal wise
	while (!is_ali.eof()) 
	{

		// read the diagonal
		Diagonal diagonal;
		is_ali >> diagonal;
		is_ali.ignore();		// skip colon

		// for a new diagonal, position yourself at the first residue 
		// on the diagonal
		if (diagonal < 0) 
		{
			row = -diagonal + row_from;
			col = 0;
		} else {
			row = 0;
			col = diagonal + col_from;
		}

		while (is_ali.peek() != ';' && !is_ali.eof()) 
		{
			Position d = 0;   
			is_ali >> d;
			// write a gap
			if (d < 0) 
			{
				row -= d;
				col -= d;
			} 
			else 
			{
				// emit dots
				while (d > 0) 
				{
					if (reverse)
						dest->addPair( col++, row++, 0 );       
					else 
						dest->addPair( row++, col++, 0 );       
					d--;
				}
			}	
		} 
		is_ali.ignore();		// skip semicolon
	}   

	return;	
}

//--------------------------------------------------------------------------------------------------------------------------------
std::ostream & operator<< (std::ostream & output, const AlignmentFormatDiagonals & src) 
{
	output 
		<< src.mRowFrom << "\t" << src.mRowTo << "\t" << "\t" 
		<< src.mColFrom << "\t" << src.mColTo << "\t" << mAlignment << "\t"
	return output;
}

std::istream & operator>> (std::istream & input, AlignmentFormatDiagonals & dest) 
{
	input >> dest.mRowFrom >> dest.mRowTo >> dest.mColFrom >> dest.mColTo >> dest.mAlignment;
	return input;
}

//----------------------------------------------------------------------------------------------------------------------
/** writes a compressed alignment in diagonal format to string out.

    Although other alignment types can be written in this format, it makes most sense
    to use it for ImplAlignmentMatrixDiagonal, because then most compression is achieved.
    In ImplAlignmentMatrixDiagonal the aligned pairs are sorted first by diagonal and then
    by row.

    @param src		alignment that is to be written
    @param	 diagonal_from	only look at diagonals larger than diagonal_from
    @param	 diagonal_to	only look at diagonals smaller than diagonal_to

 */

} // namespace alignlib
