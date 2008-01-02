/*
  alignlib - a library for aligning protein sequences

  $Id: ImplAlignandum.cpp,v 1.2 2004/01/07 14:35:33 aheger Exp $

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
#include <iomanip>
#include <cassert>
#include "alignlib.h"
#include "AlignlibDebug.h"
#include "AlignException.h"
#include "ImplAlignandum.h"
#include "Translator.h"
#include "HelpersTranslator.h"

using namespace std;

namespace alignlib 
{

//--------------------------------------------------------------------------------------
ImplAlignandum::ImplAlignandum( const Translator * translator ) : 
	mFrom(NO_POS), mTo(NO_POS), mLength(0), mIsPrepared(false), 
	mTranslator(translator) 
{
	if (mTranslator == NULL) 
		mTranslator = getDefaultTranslator();
}    

//--------------------------------------------------------------------------------------
ImplAlignandum::~ImplAlignandum () 
{
	debug_func_cerr(5);
}

//--------------------------------------------------------------------------------------
ImplAlignandum::ImplAlignandum(const ImplAlignandum & src) 
{
	debug_func_cerr(5);

	mFrom   = src.mFrom;
	mTo     = src.mTo;
	mLength = src.mLength;
	mIsPrepared = src.mIsPrepared;
	mTranslator = src.mTranslator;
}

//--------------------------------------------------------------------------------------
void ImplAlignandum::mask( Position from, Position to) 
{
	Position j;
	for (j = from; j < to; j++) 
		mask( j );
}

//--------------------------------------------------------------------------------------
Position ImplAlignandum::getLength() const 
{
	return (mTo - mFrom);
}

//--------------------------------------------------------------------------------------
void ImplAlignandum::useSegment(Position from, Position to) 
{
	assert( from <= to );

	if (from == NO_POS)
		mFrom = 0;
	else
		mFrom = from;

	if (to == NO_POS)
		mTo = mLength;
	else
		mTo = std::min( to, mLength );
}

//--------------------------------------------------------------------------------------
const Translator * ImplAlignandum::getTranslator() const 
{
	return mTranslator;
}

//--------------------------------------------------------------------------------------
Position ImplAlignandum::getFrom() const 
{
	return mFrom;
}

//--------------------------------------------------------------------------------------
Position ImplAlignandum::getTo() const 
{
	return mTo;
}

//--------------------------------------------------------------------------------------
void ImplAlignandum::setTrueLength( Position length) const 
{
	mLength = length;
}

//--------------------------------------------------------------------------------------
Position ImplAlignandum::getFullLength() const 
{
	return mLength;
}

//--------------------------------------------------------------------------------------
bool ImplAlignandum::isPrepared() const 
{
	return mIsPrepared;
}

//--------------------------------------------------------------------------------------
void ImplAlignandum::setPrepared( bool flag ) const 
{
	mIsPrepared = flag;
} 

//--------------------------------------------------------------------------------------
char ImplAlignandum::asChar( Position pos ) const 
{
	return mTranslator->decode( asResidue( pos ));
}


//--------------------------------------------------------------------------------------
// use faster implementations in subclasses, if you prefer
std::string ImplAlignandum::asString() const 
{
	std::string ret_val("");

	for (Position i = 0; i < mLength; i++) 
		ret_val += mTranslator->decode( asResidue(i) );

	return ret_val;
}

void ImplAlignandum::save( std::ostream & output ) const
{
	__save( output, MNNoType );
}

//--------------------------------------------------------------------------------------
void ImplAlignandum::__save( std::ostream & output, MagicNumberType type ) const 
{
	if (type == MNNoType )
	{
		type = MNImplAlignandum;
		output.write( (char*)&type, sizeof(MagicNumberType) );
	}

	output.write( (char*)&mIsPrepared, sizeof(bool) );
	output.write( (char*)&mFrom, sizeof(Position) );
	output.write( (char*)&mTo, sizeof(Position) );
	output.write( (char*)&mLength, sizeof(Position) );
}

//--------------------------------------------------------------------------------------
void ImplAlignandum::load( std::istream & input)  
{

	input.read( (char*)&mIsPrepared, sizeof(bool) );
	input.read( (char*)&mFrom, sizeof(Position) );
	input.read( (char*)&mTo, sizeof(Position) );
	input.read( (char*)&mLength, sizeof(Position) );

	if (input.fail()) 
		throw AlignException( "incomplete Alignandum object in stream.");
	
}


} // namespace alignlib
