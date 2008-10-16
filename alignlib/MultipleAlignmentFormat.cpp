//--------------------------------------------------------------------------------
// Project LibAlign
//
// Copyright (C) 2000 Andreas Heger All rights reserved
//
// Author: Andreas Heger <heger@ebi.ac.uk>
//
// $Id: MultipleAlignmentFormat.cpp,v 1.1.1.1 2002/07/08 21:20:17 heger Exp $
//--------------------------------------------------------------------------------    


#include <iostream>
#include <iomanip>
#include <iterator>
#include <cstring>
#include <string>
#include <sstream>

#include "alignlib_types.h"
#include "alignlib_interfaces.h"
#include "AlignlibDebug.h"
#include "AlignlibException.h"
#include "MultipleAlignmentFormat.h"
#include "HelpersEncoder.h"
#include "HelpersAlignment.h"
#include "HelpersAlignatum.h"

using namespace std;

namespace alignlib 
{

//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
MultipleAlignmentFormat::MultipleAlignmentFormat() : mRepresentation("")
{
}

MultipleAlignmentFormat::MultipleAlignmentFormat( const MultipleAlignmentFormat & src) :
	mRepresentation(src.mRepresentation)
{
}

MultipleAlignmentFormat::MultipleAlignmentFormat( std::istream & input ) 
{
	load( input );
}

MultipleAlignmentFormat::MultipleAlignmentFormat( const std::string & src) 
{
	std::istringstream i(src.c_str());
	load( i );
}

MultipleAlignmentFormat::MultipleAlignmentFormat( const HMultipleAlignment & src) 
{
	fill( src );
}

MultipleAlignmentFormat::~MultipleAlignmentFormat()
{
}

void MultipleAlignmentFormat::fill( const HMultipleAlignment & src )
{
	debug_func_cerr( 5 );
}

void MultipleAlignmentFormat::copy( HMultipleAlignment & dest ) const
{
	debug_func_cerr( 5 );
	dest->clear();
}

void MultipleAlignmentFormat::load( std::istream & input)
{
	debug_func_cerr( 5 );
	input >> mRepresentation;
}

void MultipleAlignmentFormat::save( std::ostream & output) const
{
	debug_func_cerr( 5 );
	output << mRepresentation;
}

//--------------------------------------------------------------------------------------------------------------------------------
std::ostream & operator<< (std::ostream & output, const MultipleAlignmentFormat & src)
{
	src.save( output );
	return output;
}

//--------------------------------------------------------------------------------------------------------------------------------
std::istream & operator>> (std::istream & input, MultipleAlignmentFormat & dest) 
{
	dest.load( input );
	return input;
}

//--------------------------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------------------
MultipleAlignmentFormatPlain::MultipleAlignmentFormatPlain() 
: MultipleAlignmentFormat()
{
}

MultipleAlignmentFormatPlain::MultipleAlignmentFormatPlain( std::istream & input ) 
: MultipleAlignmentFormat()
{
	load( input );
}

MultipleAlignmentFormatPlain::MultipleAlignmentFormatPlain( const std::string & src) 
: MultipleAlignmentFormat()
	{
	std::istringstream i(src.c_str());
	load( i );
	}

MultipleAlignmentFormatPlain::MultipleAlignmentFormatPlain( const HMultipleAlignment & src) 
: MultipleAlignmentFormat()
{
	fill( src );
}


MultipleAlignmentFormatPlain::~MultipleAlignmentFormatPlain () 
{
}

MultipleAlignmentFormatPlain::MultipleAlignmentFormatPlain (const MultipleAlignmentFormatPlain & src ) 
: MultipleAlignmentFormat( src )
{
}

void MultipleAlignmentFormatPlain::fill( const HMultipleAlignment & src)
{
	debug_func_cerr(5);

	MultipleAlignmentFormat::fill( src );
	for (int x = 0; x < src->getNumSequences(); ++x)
	{
		mRepresentation += src->getRow(x)->getString() + '\n'; 
	}
	
}

//--------------------------------------------------------------------------------------------------------------------------------
void MultipleAlignmentFormatPlain::copy( HMultipleAlignment & dest ) const 
{
	debug_func_cerr(5);

	MultipleAlignmentFormat::copy( dest );
	assert( false );
}


} // namespace alignlib
