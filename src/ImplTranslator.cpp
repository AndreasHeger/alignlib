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
#include "alignlib.h"
#include "AlignlibDebug.h"
#include "ImplTranslator.h"

using namespace std;

namespace alignlib 
{

//---------------------------------------------------------< constructors and destructors >--------------------------------------
ImplTranslator::ImplTranslator () : 
	mAlphabet( ""), mEncodingTable(0) 
	{
	}

ImplTranslator::~ImplTranslator () 
{
	if ( mEncodingTable != NULL )
		delete [] mEncodingTable;
}

ImplTranslator::ImplTranslator (const ImplTranslator & src ) :
	mAlphabet( src.mAlphabet ),
	mEncodingTable( NULL ) 
	{
	if (src.mEncodingTable != NULL)
	{
		mEncodingTable = new Residue[ sizeof( char ) ];
		memcpy( mEncodingTable, 
				src.mEncodingTable, 
				sizeof(char) * mAlphabet.size() );		
	}
	}

ImplTranslator::ImplTranslator(  const std::string & alphabet ) : 
	mAlphabet( alphabet ), mEncodingTable( NULL ) 
	{
	mEncodingTable = new Residue[ sizeof( char) ];

	for ( int x = 0; x < sizeof(char); ++x )
		mEncodingTable[x] = CODE_MASK;

	for ( int x = 0; x < mAlphabet.size() ; ++x)
	{
		mEncodingTable[(Residue)toupper(mAlphabet[x])] = x;
		mEncodingTable[(Residue)tolower(mAlphabet[x])] = x;
	}
	}

//--------------------------------------------------------------------------------------------------------------------------------
char ImplTranslator::decode( const Residue residue) const 
{ 
	return ( (residue == CODE_GAP) ? CHAR_GAP : mAlphabet[residue] ); 
}

//--------------------------------------------------------------------------------------------------------------------------------
std::string ImplTranslator::decode( const Residue * src, int length ) const 
{
	debug_func_cerr(5);

	char * result = new char[length + 1];

	int i;
	for (i = 0; i < length; i++) 
		result[i] = (src[i] == CODE_GAP) ? CODE_GAP : mAlphabet[src[i]]; 

	result[length] = '\0'; 
	std::string s( result );
	delete [] result;
	return s;
}    

//--------------------------------------------------------------------------------------------------------------------------------
Residue ImplTranslator::encode( const char residue) const 
{ 
	return mEncodingTable[residue]; 
}  

//--------------------------------------------------------------------------------------------------------------------------------
Residue * ImplTranslator::encode( const std::string & src ) const 
{

	Residue * result = new Residue[src.size()];

	int i;
	for (i = 0; i < src.size(); i++)
		result[i] = mEncodingTable[src[i]];

	return result;
}

//--------------------------------------------------------------------------------------------------------------------------------
bool ImplTranslator::isValidChar( const char query ) const 
{
	// todo: use std::string functions
	for ( int x = 0; x < mAlphabet.size(); ++x)
	{
		if ( mAlphabet[x] == query )
			return true;
	}
	return false;
}    

//--------------------------------------------------------------------------------------------------------------------------------
int ImplTranslator::getAlphabetSize() const
{
	return mAlphabet.size();
}


//--------------------------------------------------------------------------------------------------------------------------------
Residue ImplTranslator::getMaskCode() const 
{
	return CODE_MASK;
}

//--------------------------------------------------------------------------------------------------------------------------------
Residue ImplTranslator::getGapCode() const 
{
	return CODE_GAP;
}

//--------------------------------------------------------------------------------------------------------------------------------
char ImplTranslator::getMaskChar() const 
{
	return CHAR_MASK;
}

//--------------------------------------------------------------------------------------------------------------------------------
char ImplTranslator::getGapChar() const 
{
	return CHAR_GAP;
}


//--------------------------------------------------------------------------------------------------------------------------------

} // namespace alignlib
