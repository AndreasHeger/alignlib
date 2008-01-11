/*
  alignlib - a library for aligning protein sequences

  $Id: HelpersTranslator.cpp,v 1.2 2004/01/07 14:35:33 aheger Exp $

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
#include "alignlib.h"
#include "alignlib_fwd.h"
#include "alignlib_default.h"
#include "AlignlibDebug.h"
#include "AlignException.h"
#include "Translator.h"
#include "ImplTranslator.h"
#include "HelpersTranslator.h"

namespace alignlib 
{

//-------------------------------------------------------------------------------

const HTranslator makeTranslator( const AlphabetType & alphabet_type )
{
	debug_func_cerr( 5 );
	ImplTranslator * t;
	switch (alphabet_type) 
	{
	case Protein20: 
		t = new ImplTranslator( Protein20, "ACDEFGHIKLMNPQRSTVWY", "-.", "X" ); 
		break;
	case Protein23:
		t = new ImplTranslator( Protein23, "ABCDEFGHIKLMNPQRSTVWXYZ", "-.", "X" );
		break;
	case DNA4: 
		t = new ImplTranslator( DNA4, "ACGT", "-.", "N" );
		break;
	default:
		throw AlignException( "unknown alphabet" );
	}
	return HTranslator( t );
}

/** get a built-in translator object 
 * */
const HTranslator getTranslator( const AlphabetType & alphabet_type )
{
	// The static variables are initialized the first time this function is
	// called and then retain their values.
	
	// 20-letter alphabet plus X
	static const HTranslator translator_protein_20(makeTranslator( Protein20)); 

	// encoding table compatible with BLOSUM and PAML matrices
	static const HTranslator translator_protein_23(makeTranslator( Protein23));

	// 4-letter DNA alphabet
	static const HTranslator translator_dna_4(makeTranslator( DNA4));
	
	debug_func_cerr( 5 );
	switch (alphabet_type) 
	{
	case Protein20: 
		return translator_protein_20; break;		
	case Protein23:
		return translator_protein_23; break;
	case DNA4: 
		return translator_dna_4; break;
	}
	throw AlignException( "unknown alphabet" );
}


/** load a translator object from stream
 */
const HTranslator loadTranslator( std::istream & input )
{
	// read Alignandum type
	AlphabetType alphabet_type;

	if (input.eof()) 
		throw AlignException("HelpersTranslator.cpp: incomplete translator.");

	input.read( (char*)&alphabet_type, sizeof(AlphabetType) );

	if (input.eof()) 
		throw AlignException("HelpersTranslator.cpp: incomplete translator - could not read alphabet type.");

	HTranslator result;

	switch (alphabet_type)
	{
	case User : 
	{
		// read user alphabet
		size_t size;
		input.read( (char *)&size, sizeof( size_t ));
		char * alphabet = new char[size];
		input.read( alphabet, sizeof(char) * size);

		input.read( (char *)&size, sizeof( size_t ));
		char * gap_chars = new char[size];
		input.read( gap_chars, sizeof(char) * size);

		input.read( (char *)&size, sizeof( size_t ));
		char * mask_chars = new char[size];
		input.read( mask_chars, sizeof(char) * size);

		if (input.eof())
			throw AlignException( "HelpersTranslator.cpp: incomplete translator ");

		result = HTranslator( new ImplTranslator( alphabet_type, alphabet, gap_chars, mask_chars ) );

		delete [] alphabet;
		delete [] gap_chars;
		delete [] mask_chars;

		break;
	}
	case Protein20 :
		result = getTranslator( Protein20 );
		break;
	case Protein23:
		result = getTranslator( Protein23 );
		break;
	case DNA4:
		result = getTranslator( DNA4 );
		break;
	default:
		throw AlignException( "HelpersTranslator: unknown object found in stream" );
	}	
	return result;
}	

} // namespace alignlib
