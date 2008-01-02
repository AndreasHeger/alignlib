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
#include "alignlib.h"
#include "AlignlibDebug.h"

#include "HelpersAlignandum.h"
#include "HelpersSequence.h"
#include "ImplSequence.h" 
#include "AlignException.h"
#include "Translator.h"
#include "HelpersTranslator.h"

using namespace std;

namespace alignlib 
{

//---------------------------------< implementation of factory functions >--------------

//--------------------------------------------------------------------------------------
ImplSequence::ImplSequence( const Translator * translator ) :
	ImplAlignandum( translator ),
	mSequence(NULL) 
{
}

//--------------------------------------------------------------------------------------
ImplSequence::ImplSequence( const std::string & src, 
		const Translator * translator  ) : 
	ImplAlignandum( translator ), mSequence(NULL) 
	{
	Position length = src.size();

	setTrueLength( length );
	useSegment();
	mSequence = translator->encode( src );
	setPrepared(true );
	}

//--------------------------------------------------------------------------------------
ImplSequence::ImplSequence( const ImplSequence & src ) : 
	ImplAlignandum( src ), mSequence(NULL)
{
	debug_func_cerr(5);

	if (mSequence != NULL) delete [] mSequence;
	//!! make exception safe
	mSequence = new Residue[src.getFullLength()];
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
Alignandum * ImplSequence::getClone() const 
{
	return new ImplSequence( *this );
}


//--------------------------------------------------------------------------------------
Residue ImplSequence::asResidue(Position n) const 
{ 
	return mSequence[n]; 
}

//--------------------------------------------------------------------------------------
const AlignandumDataSequence & ImplSequence::getData() const 
{
	mData.mSequencePointer = mSequence;
	return mData;
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
void ImplSequence::mask( Position x) 
{
	mSequence[ x ] = mTranslator->getMaskCode();
}


//--------------------------------------------------------------------------------------
void ImplSequence::shuffle( unsigned int num_iterations, Position window_size ) 
{

	if (window_size == 0)
		window_size = getLength();

	Position first_from = getFrom();

	for (unsigned x = 0; x < num_iterations; x++) 
	{ 

		Position i,j;
		Position to = getTo();

		while (to > first_from ) 
		{
			Position from = to - window_size;

			if (from < 1) 
			{
				from = 1;
				window_size = to;
			}

			for (i = to; i > from; i--) 
			{
				j = to - getRandomPosition(window_size) - 1;
				Residue x = mSequence[j];
				mSequence[j] = mSequence[i];
				mSequence[i] = x;
			}

			to -= window_size;
		}
	}
}


//--------------------------------------------------------------------------------------
void ImplSequence::write( std::ostream & output ) const 
{
	std::string s = mTranslator->decode( mSequence, getFullLength() );
	output << s;
}

//--------------------------------------------------------------------------------------
void ImplSequence::__save( std::ostream & output, MagicNumberType type ) const 
{
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
	ImplAlignandum::load( input );

	if (mSequence != NULL) 
		delete [] mSequence;
	
	mSequence = new Residue[ getFullLength() ] ;
	input.read( (char*)mSequence, sizeof(Residue) * getFullLength() );
	
	if (input.fail()) 
		throw AlignException( "incomplete sequence in stream.");
	
}



} // namespace alignlib
