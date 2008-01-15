/*
  alignlib - a library for aligning protein sequences

  $Id: ImplTranslator.cpp,v 1.4 2004/09/16 16:02:38 aheger Exp $

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


#include <string.h>
#include <iostream>
#include "alignlib_fwd.h"
#include "alignlib_interfaces.h"
#include "AlignlibDebug.h"
#include "AlignException.h"
#include "ImplTranslator.h"

using namespace std;

namespace alignlib 
{

//---------------------------------------------------------< constructors and destructors >--------------------------------------
ImplTranslator::ImplTranslator () : 
	mAlphabetType( User ), mAlphabet( ""), mGapChars( "" ), mMaskChars(""), 
	mTableSize(0), mEncodingTable(0), mDecodingTable(0), mAlphabetSize( 0 ) 
	{
	debug_func_cerr( 5 );
	}

//--------------------------------------------------------------------------------------------------------------------------------
ImplTranslator::~ImplTranslator () 
{
	debug_func_cerr( 5 );
	if ( mEncodingTable != NULL )
		delete [] mEncodingTable;
	
	if ( mDecodingTable != NULL )
		delete [] mDecodingTable;
}

//--------------------------------------------------------------------------------------------------------------------------------
ImplTranslator::ImplTranslator (const ImplTranslator & src ) :
	mAlphabet( src.mAlphabet ),
	mTableSize( src.mTableSize ),
	mEncodingTable( NULL ),
	mDecodingTable( NULL ),
	mAlphabetSize( src.mAlphabetSize )
{
	debug_func_cerr( 5 );
	
	if (src.mEncodingTable != NULL)
	{
		mEncodingTable = new Residue[ mTableSize ];
		memcpy( mEncodingTable, 
				src.mEncodingTable, 
				sizeof(char) * mTableSize );
		mDecodingTable = new char[ mTableSize ];
		memcpy( mDecodingTable, 
				src.mDecodingTable, 
				sizeof(char) * mTableSize );
		
		
	}
	// this might confuse built-in versus non-builtin objects
	// TODO: check built-in versus non-built-in behaviour
	mAlphabetType = src.mAlphabetType;
}

//--------------------------------------------------------------------------------------------------------------------------------
ImplTranslator::ImplTranslator ( const AlphabetType & alphabet_type,
								 const std::string & alphabet,
								 const std::string & gap_chars,
								 const std::string & mask_chars ) :
									 mAlphabetType( alphabet_type ), 
									 mAlphabet( alphabet ), 
									 mGapChars( gap_chars ),
									 mMaskChars( mask_chars ),
									 mEncodingTable( NULL ),
									 mAlphabetSize( 0 )
{
	
	debug_func_cerr( 5 );
	
	// assertions to check for empty input
	if (mGapChars.size() == 0)
		throw AlignException( "ImplTranslator.cpp: no gap characters specified.");
	
	if (mMaskChars.size() == 0)
		throw AlignException( "ImplTranslator.cpp: no mask characters specified.");
	
	if (mAlphabet.size() == 0 )
		throw AlignException( "ImplTranslator.cpp: alphabet is empty.");
	
	// build encoding and decoding table
	mTableSize = std::numeric_limits<char>::max();
		
	mEncodingTable = new Residue[ mTableSize + 1 ];
	mDecodingTable = new char[ mTableSize + 1];

	for ( Residue x = 0; x <= mTableSize; ++x )
	{
		mEncodingTable[x] = mTableSize;
		mDecodingTable[x] = mMaskChars[0];
	}
	 mAlphabetSize = 0;
	
	for ( Residue x = 0; x < mAlphabet.size() ; ++x)
	{
		mEncodingTable[(unsigned int)toupper(mAlphabet[x])] = mAlphabetSize;
		mEncodingTable[(unsigned int)tolower(mAlphabet[x])] = mAlphabetSize;
		mDecodingTable[mAlphabetSize] = mAlphabet[x];
		++mAlphabetSize;
	}
	
	Residue mask_code = mAlphabetSize;
	char mask_char = mMaskChars[0]; 
	
	// masking characters can appear in the alphabet (to ensure they use a specific index)
	for ( Residue x = 0; x < mMaskChars.size() ; ++x)
	{
		if (mEncodingTable[mMaskChars[x]] == mTableSize )
		{
			mEncodingTable[(unsigned int)toupper(mMaskChars[x])] = mask_code;
			mEncodingTable[(unsigned int)tolower(mMaskChars[x])] = mask_code;
			mDecodingTable[mAlphabetSize] = mask_char;
			++mAlphabetSize;
		}
	}
		
	// set all unknown characters to the masking character
	for ( Residue x = 0; x <= mTableSize; ++x )
		if (mEncodingTable[x] == mTableSize)
			mEncodingTable[x] = mask_code;
	
	// map gap characters to maximum index
	for ( Residue x = 0; x < mGapChars.size(); ++x)
		mEncodingTable[(unsigned int)mGapChars[x]] = mTableSize;
	mDecodingTable[mTableSize] = mGapChars[0];
}

//--------------------------------------------------------------------------------------------------------------------------------
char ImplTranslator::operator[]( const Residue & residue ) const
{
	return mDecodingTable[residue]; 
}

//--------------------------------------------------------------------------------------------------------------------------------
Residue ImplTranslator::operator[]( const char & c ) const
{
	return mEncodingTable[(unsigned int)c]; 	
}

//--------------------------------------------------------------------------------------------------------------------------------
char ImplTranslator::decode( const Residue residue) const 
{ 
	return mDecodingTable[residue]; 
}

//--------------------------------------------------------------------------------------------------------------------------------
std::string ImplTranslator::decode( const Residue * src, int length ) const 
{
	debug_func_cerr(5);

	char * result = new char[length + 1];

	int i;
	for (i = 0; i < length; i++) 
		result[i] = mDecodingTable[src[i]];

	result[length] = '\0'; 
	std::string s( result );
	delete [] result;
	return s;
}    

//--------------------------------------------------------------------------------------------------------------------------------
Residue ImplTranslator::encode( const char residue) const 
{ 
	return mEncodingTable[(unsigned int)residue]; 
}  

//--------------------------------------------------------------------------------------------------------------------------------
HResidueVector ImplTranslator::encode( const std::string & src ) const 
{

	HResidueVector result(new ResidueVector( src.size() ));
	for ( int i = 0; i < src.size(); i++)
		(*result)[i] = mEncodingTable[src[i]];

	return result;
}

//--------------------------------------------------------------------------------------------------------------------------------
bool ImplTranslator::isValidChar( const char query ) const 
{
	return ( mAlphabet.find( query ) != std::string::npos ||
			mMaskChars.find( query ) != std::string::npos );
}    

//--------------------------------------------------------------------------------------------------------------------------------
int ImplTranslator::getAlphabetSize() const
{
	return mAlphabetSize ;
}

//--------------------------------------------------------------------------------------------------------------------------------
std::string ImplTranslator::getAlphabet() const
{
	return mAlphabet;
}

//--------------------------------------------------------------------------------------------------------------------------------
AlphabetType ImplTranslator::getAlphabetType() const
{
	return mAlphabetType;
}

//--------------------------------------------------------------------------------------------------------------------------------
Residue ImplTranslator::getMaskCode() const 
{
	return encode( mMaskChars[0] );
}

//--------------------------------------------------------------------------------------------------------------------------------
Residue ImplTranslator::getGapCode() const 
{
	return encode( mGapChars[0] );
}

//--------------------------------------------------------------------------------------------------------------------------------
std::string ImplTranslator::getMaskChars() const 
{
	return mMaskChars;
}

//--------------------------------------------------------------------------------------------------------------------------------
std::string ImplTranslator::getGapChars() const 
{
	return mGapChars;
}

//--------------------------------------------------------------------------------------------------------------------------------
char ImplTranslator::getMaskChar() const 
{
	return mMaskChars[0];
}

//--------------------------------------------------------------------------------------------------------------------------------
char ImplTranslator::getGapChar() const 
{
	return mGapChars[0];
}

//--------------------------------------------------------------------------------------
void ImplTranslator::write( std::ostream & output ) const 
{
	for ( Residue x = 0; x < mAlphabet.size(); ++x )
		output << (int)x << '\t' << mAlphabet[x] << '\t' << (int)encode(mAlphabet[x]) << '\t' << decode(encode(mAlphabet[x]))<< std::endl;

	output << getGapChar() << '\t' << (int)getGapCode() << std::endl;
	output << getMaskChar() << '\t' << (int)getMaskCode() << std::endl;
}
	
//--------------------------------------------------------------------------------------
void ImplTranslator::save( std::ostream & output ) const 
{
	debug_func_cerr( 5 );
	output.write( (char *)&mAlphabetType, sizeof( AlphabetType) );
	
	if ( mAlphabetType == User )
	{	
		output.write( (char *)mAlphabet.size(), sizeof( size_t ) );
		output.write( (char *)mAlphabet.c_str(), mAlphabet.size() * sizeof( char ) );
		output.write( (char *)mGapChars.size(), sizeof( size_t ) );		
		output.write( (char *)mGapChars.c_str(), mGapChars.size() * sizeof( char ) );
		output.write( (char *)mMaskChars.size(), sizeof( size_t ) );		
		output.write( (char *)mMaskChars.c_str(), mMaskChars.size() * sizeof( char ) );
	}
}

//--------------------------------------------------------------------------------------
// This will map alphabet and mask characters, but not gap characters.
//
HResidueVector ImplTranslator::map( const HTranslator & other ) const
{
	debug_func_cerr( 5 );
	HResidueVector map_other2this( new ResidueVector( other->getAlphabetSize(), getMaskCode()) );

	for ( Residue x = 0; x < other->getAlphabetSize(); ++x)
		(*map_other2this)[x] = encode( other->decode( x ) );

	return map_other2this;
}



//--------------------------------------------------------------------------------------------------------------------------------

} // namespace alignlib
