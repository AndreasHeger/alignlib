/*
  alignlib - a library for aligning protein sequences

  $Id: test_MultipleAlignment.cpp,v 1.6 2004/06/02 12:11:38 aheger Exp $

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

/** Test the Translator object
 */

#if HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>
#include <fstream>
#include <sstream>
#include <cassert>

#include "Translator.h"
#include "HelpersTranslator.h"

using namespace std;
using namespace alignlib;

void testTranslator( const HTranslator & translator, 
		const std::string & alphabet,
		const std::string & gap_chars,
		const std::string & mask_chars,
		int alphabet_size )
{
	std::cout << *translator << std::endl;
	
	// check alphabet size (number of chars + 1 for mask value)
	assert( alphabet_size == translator->getAlphabetSize() );
	
	for (int x = 0; x < alphabet.size(); ++x)
	{
		assert( x == translator->encode( alphabet[x] ) );
		// check if upper/lower-case is ignored
		assert( translator->encode( toupper(alphabet[x]) ) == translator->encode( tolower( alphabet[x] ) ) );
	}

	// check mask codes
	for (int x = 0; x < mask_chars.size(); ++ x)
		assert( translator->encode(mask_chars[x]) == translator->getMaskCode() );

	// check gap codes
	for (int x = 0; x < gap_chars.size(); ++ x)
		assert( translator->encode(gap_chars[x]) == translator->getGapCode() );

	// check saving/loading from stream.
	{
		ofstream file("test_Translator.tmp", ios::binary);
		translator->save( file );
		translator->save( file );		
		file.close();
	}
	
	{
		ifstream file("test_Translator.tmp", ios::binary) ;
		
		HTranslator b;
		int n = 0;
		while ( file.peek() != EOF )
		{
			b = loadTranslator( file ) ;
			if (b->getAlphabetType() != User )
				assert( translator == b );
			++n; 
		}
		assert( n == 2 );
	}	
	
	{
		HTranslator n( getTranslator( Protein20 ));
		HResidueVector map_new2old( n->map( translator ) );
		for (Residue x = 0; x < map_new2old->size(); ++x)
			std::cout << (int)x << "\t" << (int)(*map_new2old)[x] << std::endl;
	}

	{
		HTranslator n( getTranslator( Protein20 ));
		HResidueVector map_new2old( translator->map( n ) );
		for (Residue x = 0; x < map_new2old->size(); ++x)
			std::cout << (int)x << "\t" << (int)(*map_new2old)[x] << std::endl;
	}
	
}


int main () 
{
	std::cout << "----- testing Protein23 ---------" << std::endl;
	{
		testTranslator( getTranslator( Protein23 ), "ABCDEFGHIKLMNPQRSTVWXYZ", "-.", "X", 23 );
	}

	std::cout << "----- testing Protein20 ---------" << std::endl;
	{
		testTranslator( getTranslator( Protein20 ), "ACDEFGHIKLMNPQRSTVWY", "-.", "X", 21 );
	}

	std::cout << "----- testing DNA4 ---------" << std::endl;
	{
		testTranslator( getTranslator( DNA4 ), "ACGT", "-.", "N", 5 );
	}
	
	
	
	
}

