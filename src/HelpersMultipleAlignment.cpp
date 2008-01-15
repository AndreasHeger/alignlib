/*
  alignlib - a library for aligning protein sequences

  $Id: HelpersMultipleAlignment.cpp,v 1.5 2004/03/19 18:23:40 aheger Exp $

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
#include <iomanip>
#include <vector>
#include <algorithm>

#include "HelpersProfile.h"
#include "ImplMultipleAlignment.h"
#include "ImplMultipleAlignmentDots.h"
#include "AlignException.h"
#include "Alignatum.h"
#include "HelpersAlignatum.h"
#include "Alignandum.h"
#include "Alignment.h"
#include "AlignmentIterator.h"
#include "HelpersMultipleAlignment.h"

#include "Regularizor.h"
#include "HelpersRegularizor.h"

#include "Translator.h"
#include "HelpersTranslator.h"

#include "Matrix.h"
#include <math.h>

using namespace std;

namespace alignlib 
{

/** factory functions */

/** create an empty multiple alignment */
HMultipleAlignment makeMultipleAlignment() 
{
	return HMultipleAlignment( new ImplMultipleAlignment() );
}

HMultipleAlignment makeMultipleAlignmentDots( bool compress_unaligned_columns,
		int max_insertion_length) 
{
	return HMultipleAlignment( new ImplMultipleAlignmentDots( compress_unaligned_columns, max_insertion_length ) );
}

//---------------------------------------------------------------------
/** extract a multiple Alignment object from a stream in FASTA format 
    uses the extraction routine from ImplSequence objects
 */
/* 
HMultipleAlignment extractMultipleAlignmentFasta( HMultipleAlignment dest, 
						   std::istream & input ) {
  dest->clear();

  while (!input.eof()) {
    Alignatum * a = makeAlignatumFasta();
    input >> *a;
    dest->add( a );
  }

  return dest;
}
 */
void fillMultipleAlignment( 
		HMultipleAlignment & ali, 
		const std::string & sequences, 
		int nsequences ) 
{

	ali->clear();

	int total_length = sequences.size();

	int length = total_length / nsequences;

	char * buffer = new char[length + 1];

	for (int i = 0; i < total_length; i+= length) 
	{
		memcpy( buffer, &sequences[i], length);
		buffer[length] = '\0';

		HAlignatum a(makeAlignatum( buffer ));

		if (a->getAlignedLength() != 0)
			ali->add( a );
	}

	delete [] buffer;

	return;
}

/*
HMultipleAlignment fillMultipleAlignment( HMultipleAlignment ali, const char * filename ) 
{

    ali->clear();

    ifstream fin( filename);  
    if (!fin)       
	throw AlignException("Could not open file in fillMultipleAlignment");

    while (!fin.eof()) {                          // while (fin) does not work!!       
	// read in lines at a time:
	std::string s;
	fin >> s;
	Alignatum * a = makeAlignatumFromString( s );
	if (a->getAlignedLength() != 0)
	    ali->add( a );
	else {
	    delete a;
	    break;
	}
    }

    fin.close();  

    return ali;

}

HMultipleAlignment fillMultipleAlignment( HMultipleAlignment ali, 
					   const char * filename, 
					   const Alignatum * alignatum_template) {

    ali->clear();

    ifstream fin( filename);  
    if (!fin)       
	throw AlignException("Could not open file in fillMultipleAlignment");

    while (!fin.eof()) {                          // while (fin) does not work!!       
	// get an empty copy Alignatum
	Alignatum * a = alignatum_template->getNew();

	// read from stream
	fin >> *a;
	if (a->getAlignedLength() != 0) 
	{
	    ali->add(a);
	} else {
	    // only read up to first mistake
	    delete a;
	    break;
	}
    }

    fin.close();  

    return ali;

}
 */

//----------------------------------------------------------------------
/** split a multiple alignment in two groups 
 */
void copyMultipleAlignment( 
		HMultipleAlignment & dest, 
		const HMultipleAlignment & src,
		unsigned int start_row,
		unsigned int end_row ) 
{

	unsigned int width = src->getWidth();

	if (end_row > width || end_row == 0) 
		end_row = width;

	dest->clear();

	unsigned int row;

	for (row = start_row; row < end_row; row++) 
		dest->add( src->getRow(row)->getClone() );
			

	return;
}




} // namespace alignlib

