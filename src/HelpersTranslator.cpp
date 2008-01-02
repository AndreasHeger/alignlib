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
#include "AlignlibDebug.h"
#include "AlignException.h"
#include "Translator.h"
#include "ImplTranslator.h"
#include "HelpersTranslator.h"

namespace alignlib 
{

//-------------------------------------------------------------------------------

// the built-in translator objects

/** 20-letter alphabet plus X */
static const ImplTranslator translator_protein_20 = ImplTranslator( Protein20, "ACDEFGHIKLMNPQRSTVWY", "-.", "X" );

/** encoding table compatible with BLOSUM and PAML matrices */
static const ImplTranslator translator_protein_23 = ImplTranslator( Protein23, "ABCDEFGHIKLMNPQRSTVWXYZ", "-.", "X" );

/** 4-letter DNA alphabet */
static const ImplTranslator translator_dna_4 = ImplTranslator( DNA4, "ACGT", "-.", "N" );

const Translator * DEFAULT_TRANSLATOR = &translator_protein_23;

const Translator * getTranslator( const AlphabetType & alphabet_type )
{
	switch (alphabet_type) 
	{
	case Protein20: return &translator_protein_20; break;
	case Protein23: return &translator_protein_23; break;
	case DNA4: return & translator_dna_4; break;
	}

	throw AlignException( "unknown alphabet" );
}

/** gets the default Translator object */ 
const Translator * getDefaultTranslator() 
{
	debug_func_cerr( 5 );	
	return DEFAULT_TRANSLATOR;
}

/** sets the default Translator object 
 * Only supply built-in objects, otherwise
 * memory leaks will occur.
 * */
void setDefaultTranslator( const Translator * translator ) 
{	
	debug_func_cerr( 5 );
	DEFAULT_TRANSLATOR = translator;
}

/** load a translator object from stream
 * returns NULL on EOF
 */

const Translator * loadTranslator( std::istream & input )
{
	// read Alignandum type
	AlphabetType alphabet_type;

	if (input.eof()) return NULL;

	input.read( (char*)&alphabet_type, sizeof(AlphabetType) );

	if (input.eof()) return NULL;

	const Translator * result = NULL;

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
			throw AlignException( "incomplete translator ");
					
		result = new ImplTranslator( alphabet_type, alphabet, gap_chars, mask_chars );
		
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
		throw AlignException( "unknown object found in stream" );
	}	
	return result;
}	

} // namespace alignlib
