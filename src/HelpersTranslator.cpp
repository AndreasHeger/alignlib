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

/** various alphabets */

/** 20-letter alphabet plus X */
static std::string alphabet_protein_21 = "ACDEFGHIKLMNPQRSTVWYX";

/** encoding table ompatible with BLOSUM and PAML matrices */
static std::string alphabet_protein_23 = "ABCDEFGHIKLMNPQRSTVWYXZ";	

/** 5-letter DNA alphabet */
static std::string alphabet_dna_5 = "ACGTN";	

// the built-in translator objects 
static ImplTranslator translator_protein_21 = ImplTranslator( alphabet_protein_21 );
static ImplTranslator translator_protein_23 = ImplTranslator( alphabet_protein_23 );
static ImplTranslator translator_dna_5 = ImplTranslator( alphabet_dna_5 );

const Translator * DEFAULT_TRANSLATOR = & translator_protein_23;

const Translator * getTranslator( const AlphabetType & alphabet_type )
{
	switch (alphabet_type) 
	{
	case Protein21: return &translator_protein_21; break;
	case Protein23: return &translator_protein_23; break;
	case DNA5: return & translator_dna_5; break;
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
void setDefaultTranslator( Translator * translator ) 
{	
	debug_func_cerr( 5 );
	DEFAULT_TRANSLATOR = translator;
}

} // namespace alignlib
