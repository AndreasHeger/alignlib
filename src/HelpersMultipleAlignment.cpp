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
#include "Alignata.h"
#include "AlignataIterator.h"
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
MultipleAlignment * makeMultipleAlignment() 
{
	return new ImplMultipleAlignment();
}

MultipleAlignment * makeMultipleAlignmentDots( bool compress_unaligned_columns,
		int max_insertion_length) 
{
	return new ImplMultipleAlignmentDots( compress_unaligned_columns, max_insertion_length );
}

//---------------------------------------------------------------------
/** extract a multiple Alignment object from a stream in FASTA format 
    uses the extraction routine from ImplSequence objects
 */
/* 
MultipleAlignment * extractMultipleAlignmentFasta( MultipleAlignment * dest, 
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
MultipleAlignment * fillMultipleAlignment( MultipleAlignment * ali, const char * sequences, int nsequences ) 
{

	ali->clear();

	int total_length = strlen(sequences);

	int length = total_length / nsequences;

	char * buffer = new char[length + 1];

	for (int i = 0; i < total_length; i+= length) 
	{

		memcpy( buffer, &sequences[i], length);
		buffer[length] = '\0';

		Alignatum * a = makeAlignatumFromString( buffer );

		if (a->getAlignedLength() != 0)
			ali->add( a );
		else 
			delete a;
	}

	delete [] buffer;

	return ali;

}

/*
MultipleAlignment * fillMultipleAlignment( MultipleAlignment * ali, const char * filename ) 
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

MultipleAlignment * fillMultipleAlignment( MultipleAlignment * ali, 
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

//------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------
std::string calculateConservation( const MultipleAlignment * mali, Frequency min_frequency) 
{

	Position row, col, length;

	Regularizor * regularizor = makeNoRegularizor();

	alignlib::Alignandum * profile = makeProfile( mali, NULL, regularizor );

	profile->prepare();

	const AlignandumDataProfile & data = (const AlignandumDataProfile &)profile->getData();
	const FrequencyColumn * frequencies = data.mFrequenciesPointer;

	length = mali->getLength();
	const Translator * translator = getDefaultTranslator();

	char * buffer = new char[length + 1];

	for (col = 1; col <= length; col++) {
		Frequency max_frequency = 0;
		Frequency f;
		Residue max_residue = translator->getGapCode();

		for (row = 0; row < PROFILEWIDTH; row++) {
			if ( (f = frequencies[col][row]) > max_frequency && f >= min_frequency) {
				max_frequency = f;
				max_residue = row;
			}
		}
		buffer[col-1] = translator->decode( max_residue );
	}

	buffer[length] = '\0';

	std::string seq(buffer);
	delete [] buffer;
	delete regularizor;
	delete profile;

	return seq;
}

//------------------------------------------------------------------------------------------------------
CountsMatrix * makeCountsByCategory( const MultipleAlignment * mali, 
		const unsigned int * map_residue2category ) {

	Position col, length;

	// build profile. Counts are calculated automatically
	Regularizor * regularizor = makeNoRegularizor();
	alignlib::Alignandum * profile = makeProfile( mali, NULL, regularizor );

	const AlignandumDataProfile & data   = (const AlignandumDataProfile &)profile->getData();
	const CountColumn    * counts = data.mCountsPointer;

	length = mali->getLength();

	// deterimine number of categories
	unsigned int num_categories;

	if (map_residue2category == NULL)
		num_categories = PROFILEWIDTH;
	else {
		num_categories = 0;
		for (unsigned int i = 0; i < PROFILEWIDTH; i++) 
			if (num_categories < map_residue2category[i]) 
				num_categories = map_residue2category[i];
	}
	num_categories++;

	// allocate and initialize result structure
	CountsMatrix * result = new CountsMatrix(length+1, num_categories);

	// go through counts and map counts to classes. Iterate
	// row-wise, so that mapping has to be done only once.
	for (unsigned int row = 0; row < PROFILEWIDTH; row++) {

		unsigned int category;
		if (map_residue2category == NULL)
			category = row;
		else
			category = map_residue2category[row];

		for (col = 1; col <= length; col++)
			(*result)[col][category] += (unsigned int)counts[col][row];
	}

	delete regularizor;
	delete profile;

	return result;
}

/** make a map from residues to categories. The following order has been suggested by Hannes for
    surface area calculations:
    'G': 0, 'P': 0
    'K': 1, 'R': 1,
    'D': 2, 'E': 2,
    'H': 3, 'F': 3, 'W':3, 'Y': 3, 'C': 3,
    'N': 4, 'Q': 4, 'S':4, 'T': 4,
    'A': 5, 'I': 5, 'L': 5,'M': 5, 'V': 5,
 */

const unsigned int MapResidue2CategorySurface[PROFILEWIDTH] = { 
		5, 3, 2, 2, 3,     /* A */
		0, 3, 5, 1, 5,     /* G */  
		5, 4, 0, 4, 1,     /* M */
		4, 4, 5, 3, 3,     /* S */
};

const unsigned int MapResidue2CategoryAll[PROFILEWIDTH] = {
		0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19 
};

const unsigned int * getMapResidue2CategorySurface() {
	return MapResidue2CategorySurface;
}

const unsigned int * getMapResidue2CategoryAll() {
	return MapResidue2CategoryAll;
}

/** return a vector of entropies calculated for a CountsMatrix
 */
VectorDouble * makeEntropyVector( const CountsMatrix * src) {

	unsigned int length      = src->getNumRows();
	unsigned int categories  = src->getNumCols();

	VectorDouble * result = new VectorDouble(length,0);

	for (unsigned int l = 0; l < length; l++) {
		double total = 0;
		for (unsigned int c = 0; c < categories; c++) {
			total += (*src)[l][c];
		}
		double e = 0;
		unsigned int counts;
		for (unsigned int c = 0; c < categories; c++) {
			if ( (counts = (*src)[l][c]) > 0) {
				double p = (double)counts / (double)total;
				e -= p * log(p);
			}
		}	
		(*result)[l] = e;
	}

	return result;
}

//----------------------------------------------------------------------
/** split a multiple alignment in two groups 
 */
MultipleAlignment * copyMultipleAlignment( MultipleAlignment * dest, 
		const MultipleAlignment * src,
		unsigned int start_row,
		unsigned int end_row ) {

	unsigned int width = src->getWidth();

	if (end_row > width || end_row == 0) 
		end_row = width;

	dest->clear();

	unsigned int row;

	for (row = start_row; row < end_row; row++) 
		dest->add( src->getRow(row)->getClone() );

	return dest;
}




} // namespace alignlib

