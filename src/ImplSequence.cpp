/*
  alignlib - a library for aligning protein sequences

  $Id: ImplSequence.cpp,v 1.2 2004/01/07 14:35:36 aheger Exp $

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
#include <string>
#include "alignlib_fwd.h"
#include "alignlib_interfaces.h"
#include "AlignlibDebug.h"

#include "HelpersAlignandum.h"
#include "HelpersSequence.h"
#include "ImplSequence.h" 
#include "AlignException.h"
#include "Encoder.h"
#include "HelpersEncoder.h"

using namespace std;

namespace alignlib 
{

//---------------------------------< implementation of factory functions >--------------

//----------------------------------------------------------------------------------
/** create a sequence from a NULL-terminated string */
HAlignandum makeSequence( const char * sequence, 
		const HEncoder & translator ) 
		{
		return makeSequence( std::string(sequence), translator );
		}

HAlignandum makeSequence( const char * sequence )
		{
		return makeSequence( std::string(sequence) );
		}

//----------------------------------------------------------------------------------
/** create a sequence from a string */
HAlignandum makeSequence( 
		const std::string & sequence,
		const HEncoder & translator ) 
		{
		return HAlignandum( new ImplSequence( sequence, translator ) );
		}

HAlignandum makeSequence( const std::string & sequence )
		{
	return HAlignandum( new ImplSequence( sequence, getDefaultEncoder() ) );
		}


//--------------------------------------------------------------------------------------
ImplSequence::ImplSequence( 
		const HEncoder & translator ) :
	ImplAlignandum( translator ),
	mSequence(NULL) 
{
}

//--------------------------------------------------------------------------------------
ImplSequence::ImplSequence( 
		const std::string & src, 
		const HEncoder & translator  ) : 
	ImplAlignandum( translator ), 
	mSequence(NULL) 
	{
	Position length = src.size();

	resize( length );
	
	for (int i = 0; i < length; ++i)
		mSequence[i] = translator->encode( src[i] );
	
	setPrepared(true );
	}

//--------------------------------------------------------------------------------------
ImplSequence::ImplSequence( const ImplSequence & src ) : 
	ImplAlignandum( src ), mSequence(NULL)
{
	debug_func_cerr(5);

	resize( src.getFullLength() );
	memcpy( mSequence, src.mSequence, src.getFullLength());
}


//--------------------------------------------------------------------------------------
ImplSequence::~ImplSequence() 
{
	debug_func_cerr(5);

	if (mSequence != NULL) 
		delete [] mSequence;
}

//--------------------------------------------------------------------------------------
void ImplSequence::resize( Position length )
{
	ImplAlignandum::resize(length);
	if (mSequence != NULL) 
		delete [] mSequence;
	mSequence = new Residue[length];
	Residue gap_code = mEncoder->getGapCode();
	
	for (Position i = 0; i < length; ++i)
		mSequence[i] = gap_code;
}

//--------------------------------------------------------------------------------------
HAlignandum ImplSequence::getClone() const 
{
	return HAlignandum( new ImplSequence( *this ) );
}


//--------------------------------------------------------------------------------------
Residue ImplSequence::asResidue(Position n) const 
{ 
	return mSequence[n]; 
}

//--------------------------------------------------------------------------------------
void ImplSequence::prepare() const 
{
}

//--------------------------------------------------------------------------------------
void ImplSequence::release() const 
{
}

//--------------------------------------------------------------------------------------
void ImplSequence::mask( const Position & x) 
{
	mSequence[ x ] = mEncoder->getMaskCode();
	ImplAlignandum::mask( x );
}

//--------------------------------------------------------------------------------------
const Residue * ImplSequence::getSequence() const 
{
	return mSequence;
}

//--------------------------------------------------------------------------------------
void ImplSequence::swap( const Position & x, const Position & y )
{
	assert( x >= 0);
	assert( x < getFullLength() );
	assert( y >= 0);
	assert( y < getFullLength() );
	std::swap( mSequence[x], mSequence[y] );
}

//--------------------------------------------------------------------------------------
void ImplSequence::write( std::ostream & output ) const 
{
	std::string s = mEncoder->decode( mSequence, getFullLength() );
	output << s;
}

//--------------------------------------------------------------------------------------
void ImplSequence::__save( std::ostream & output, MagicNumberType type ) const 
{
	debug_func_cerr( 5 );
	
	if (type == MNNoType )
	{
		type = MNImplSequence;
		output.write( (char*)&type, sizeof(MagicNumberType ) );
	}

	ImplAlignandum::__save( output, type );
	
	output.write( (char*)mSequence, sizeof(Residue) * getFullLength() );
}

//--------------------------------------------------------------------------------------
void ImplSequence::load( std::istream & input)  
{
	debug_func_cerr( 5 );
	
	ImplAlignandum::load( input );

	if (mSequence != NULL) 
		delete [] mSequence;
	
	mSequence = new Residue[ getFullLength() ] ;
	input.read( (char*)mSequence, sizeof(Residue) * getFullLength() );
	
	if (input.fail()) 
		throw AlignException( "incomplete sequence in stream.");
	
}



} // namespace alignlib
