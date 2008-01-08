/*
  alignlib - a library for aligning protein sequences

  $Id: HelpersSequence.cpp,v 1.2 2004/01/07 14:35:32 aheger Exp $

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

#include <time.h>

#include <iostream>
#include <fstream>
#include <string>

#include "Matrix.h"

#include "HelpersSequence.h" 
#include "AlignException.h"

#include "Translator.h"
#include "HelpersTranslator.h" 

#include "Alignandum.h"
#include "ImplSequence.h"

#include <math.h>

using namespace std;

namespace alignlib 
{

//---------------------------------< implementation of useful functions >--------------

//----------------------------------------------------------------------------------
/** create a sequence from a NULL-terminated string */
HAlignandum makeSequence( const char * sequence, 
		const HTranslator & translator ) 
		{
		return HAlignandum( new ImplSequence( std::string(sequence), translator ) );
		}

HAlignandum makeSequence( const char * sequence )
		{
		return makeSequence( sequence, getDefaultTranslator() );
		}

//----------------------------------------------------------------------------------
/** create a sequence from a string */
HAlignandum makeSequence( const std::string & sequence,
		const HTranslator & translator ) 
		{
		return HAlignandum( new ImplSequence( sequence, translator ) );
		}

HAlignandum makeSequence( const std::string & sequence )
		{
	return makeSequence( sequence );
		}

//----------------------------------------------------------------------------------
/** create a sequence from a stream */
/*
HAlignandum extractSequence( std::istream & input, const HTranslator & translator ) 
{
	// TODO to be implemented
	return NULL;
}
*/

//----------------------------------------------------------------------------------
/** create a sequence from a stream, put description into field description. Return Null, if unsuccessfull */
/*
HAlignandum extractSequenceFasta( std::istream & input, 
		std::string & description,
		const HTranslator & translator ) 
		{

#define MAX_CHUNK 10000

	char * buffer = new char[MAX_CHUNK];

	while ( (input.peek() != '>') && 
			!input.eof() ) {
		input.getline( buffer, MAX_CHUNK);
	}

	if (input.eof())
		return NULL;

	input.get();
	input.getline(buffer, MAX_CHUNK);

	description = buffer;
	// erase end of line
	description.erase(description.size(),1); 

	std::string sequence("");

	// build the sequence character-wise
	while ( (input.peek() != '>') && 
			!input.eof() ) {

		input.getline( buffer, MAX_CHUNK);

		for (unsigned int i = 0; i < strlen(buffer); i++) 
			if (translator->isValidChar( buffer[i] )) 
				sequence += buffer[i];
	}

	delete [] buffer;

	if (sequence.size() > 0)
		return makeSequence( sequence.c_str(), translator );
	else
		return NULL;
		}
*/
//----------------------------------------------------------------------------------
/** read a sequence from a file, given the filename */
/*
HAlignandum readSequence( const char * filename,
		const HTranslator & translator ) 
		{

	ifstream fin( filename);  

	if (!fin)       
		throw AlignException("Could not open file in ImplSequence.cpp");

	std::string sequence;

	fin >> sequence;

	fin.close();  

	return makeSequence( sequence.c_str(), translator );

		}
*/

/* this routine should be replace by something, that
     a. uses a different random generator
     b. is more efficient.
     c. is portable
 */

const double max_rand = pow(2.0,31) -1;

Residue sampleFromDistribution( const double * histogram, const int width ) 
{
	double x = random() / max_rand;
	double s = 0;
	Residue i = 0;
	for (i = 0; i < width; i++) 
	{
		s+= histogram[i];
		if (x < s) return i;
	}
	return width - 1;
}

void setRandomSeed( long seed ) 
{
	srandom(seed);
}

HAlignandum makeMutatedSequence( HAlignandum src, 
		const MutationMatrix * matrix) 
		{

	assert( matrix->getNumRows() == matrix->getNumCols() );
	// TODO: check function
	
	int width = matrix->getNumRows();
	
	// intialize random generator
	char * buffer = new char[src->getLength() + 1];

	const HTranslator & translator = getDefaultTranslator();
	unsigned int x;
	Position i;

	for (i = 0, x = 0; i < src->getLength(); i++, x++) 
	{
		Residue residue = src->asResidue(i);
		Residue new_residue = sampleFromDistribution( (*matrix)[residue],
		                                       	width );
		buffer[x] = translator->decode(new_residue);
	}

	buffer[x] = '\0';

	HAlignandum sequence = makeSequence(buffer, translator );
	delete [] buffer;

	return sequence;

		}


} // namespace alignlib















