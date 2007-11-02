/*
  alignlib - a library for aligning protein sequences

  $Id: HelpersProfile.cpp,v 1.4 2004/03/19 18:23:40 aheger Exp $

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
#include <stdio.h>

#include "alignlib.h"
#include "AlignlibDebug.h"
#include "AlignException.h"
#include "HelpersProfile.h"
#include "Weightor.h"
#include "LogOddor.h"
#include "Regularizor.h"

#include "MultipleAlignment.h"
#include "HelpersMultipleAlignment.h"

#include "Alignata.h"
#include "AlignataIterator.h"

/** default objects */
#include "Translator.h"
#include "HelpersTranslator.h"

#include "HelpersWeightor.h"
#include "HelpersRegularizor.h"
#include "HelpersLogOddor.h"

#include "ImplProfile.h"
#include "ImplSequence.h"


#ifdef WITH_DMALLOC
#include <dmalloc.h>
#endif

using namespace std;

namespace alignlib {

//---------------------------------< implementation of factory functions >--------------

typedef int TYPE_INT_COUNTS_COLUMN[PROFILEWIDTH];

//------------------------------------------------------------------------------------------
/** create empty profile */
Alignandum * makeProfile( const Regularizor * regularizor,
		const LogOddor * logoddor ) {
	if (!regularizor) 
		regularizor = getDefaultRegularizor();

	if (!logoddor)
		logoddor = getDefaultLogOddor();

	Alignandum * profile = new ImplProfile( regularizor, logoddor );

	return profile;
}

//------------------------------------------------------------------------------------------
/** create empty profile with given length */
Alignandum * makeProfile( Position length,
		const Regularizor * regularizor,
		const LogOddor * logoddor ) {
	if (!regularizor) 
		regularizor = getDefaultRegularizor();

	if (!logoddor)
		logoddor = getDefaultLogOddor();


	ImplProfile * profile = new ImplProfile( regularizor, logoddor );
	profile->setTrueLength( length );
	profile->useFullLength();
	profile->allocateCounts();
	return profile;
}



//------------------------------------------------------------------------------------------
/** create profile from a string of sequences */
Alignandum * makeProfile( const char * src, int nsequences,
		const Weightor * weightor, 
		const Regularizor * regularizor,
		const LogOddor * logoddor ) {
	if (!weightor)
		weightor = getDefaultWeightor();

	if (!regularizor) 
		regularizor = getDefaultRegularizor();

	if (!logoddor)
		logoddor = getDefaultLogOddor();

	MultipleAlignment * m = fillMultipleAlignment( makeMultipleAlignment(), src, nsequences );

	Alignandum * profile = new ImplProfile( regularizor, logoddor );
	fillProfile( profile, m );
	delete m;
	return profile;
}

//------------------------------------------------------------------------------------------
/** create a default profile from a multiple alignment */
Alignandum * makeProfile( const MultipleAlignment * mali, 
		const Weightor * weightor, 
		const Regularizor * regularizor,
		const LogOddor * logoddor ) {

	if (!weightor)
		weightor = getDefaultWeightor();

	if (!regularizor) 
		regularizor = getDefaultRegularizor();

	if (!logoddor)
		logoddor = getDefaultLogOddor();

	Alignandum * profile = new ImplProfile( regularizor, logoddor );
	fillProfile( profile, mali, weightor );
	return profile;
}

//---------------------------------------------------------------------------------------------------------------
Alignandum * fillProfile( Alignandum * dest, 
		const MultipleAlignment * src, 
		const Weightor * weightor ) 
		{
	debug_func_cerr(5);


	/* first check, that we actually do have a profile here. */
	ImplProfile * profile = dynamic_cast<ImplProfile*>(dest);

	// set up the weightor object and calculate the weights.
	if (weightor == NULL) 
		weightor = getDefaultWeightor();

	SequenceWeights * weights = weightor->calculateWeights( *src );
#ifdef DEBUG
	cout << "-------------->Weights start-----------" << endl;
	for (int i = 0; i < src->getWidth(); i++) 
		cout << i << " " << (*weights)[i] << endl;
	cout << "-------------->Weights end-------------" << endl;
#endif

	// ask profile to allocate new memory for the counts.
	// old memory is automatically released when allocating new memory, so do not worry about this here.

	Position length = src->getLength();

	profile->setTrueLength( length );
	profile->useFullLength();
	profile->allocateCounts();

	CountColumn * counts = profile->getData().mCountsPointer;

	// calculate counts
	int width = src->getWidth();

	Residue code;
	for (int nsequence = 0; nsequence < width; nsequence++) {
		const std::string & seq = (*src)[nsequence];
		for (int column = 0; column < length; column++) 
			if ( (code = getDefaultTranslator()->encode( seq[column] ) ) < PROFILEWIDTH)
				counts[column][code] += (*weights)[nsequence];
	}

	delete weights;

	profile->setPrepared( false );  

	return profile;
		}     

//------------------------------------------------------------------------------------------
/** write counts from profile in binary format to stream */
void writeProfileBinaryCounts( std::ostream & output, const Alignandum * src) {

	/* first check, that we actually do have a profile here. */
	const ImplProfile * p = dynamic_cast<const ImplProfile*>(src);

	/* If not, there is nothing to write */
	if (!p) 
		return;

	const CountColumn * counts = p->getData().mCountsPointer;

	output.write( (char*)&(counts[1]), p->getLength() * PROFILEWIDTH * sizeof( Count));

}

//------------------------------------------------------------------------------------------
/** write counts from profile in binary format to stream, store as ints with bytes bytes (not supported yet). */
void writeProfileBinaryCountsAsInt( std::ostream & output, const Alignandum * src, int bytes, float scale_factor) {


	/* first check, that we actually do have a profile here. */
	const ImplProfile * p = dynamic_cast<const ImplProfile*>(src);

	/* If not, there is nothing to write */
	if (!p) 
		return;

	const CountColumn * counts = p->getData().mCountsPointer;

	int length = p->getLength();

	int i, j;

	// allocate memory for converted profile
	TYPE_INT_COUNTS_COLUMN * int_counts = new TYPE_INT_COUNTS_COLUMN[ length + 1];

	for (i = 1; i <= length; i++)
		for (j = 0; j < PROFILEWIDTH; j++) 
			int_counts[i][j] = (int)(counts[i][j] * scale_factor);

	output.write( (char*)&(int_counts[1]), length * PROFILEWIDTH * sizeof( int ));

	delete [] int_counts;

}

//------------------------------------------------------------------------------------------
/** read counts of a profile from stream in binary format stored as integers */
Alignandum * extractProfileBinaryCountsAsInt( std::istream & input, 
		const Position max_length,
		int bytes, 
		float scale_factor,
		const Regularizor * regularizor,
		const LogOddor * logoddor ) {

	if (!regularizor) 
		regularizor = getDefaultRegularizor();

	if (!logoddor)
		logoddor = getDefaultLogOddor();

	TYPE_INT_COUNTS_COLUMN * int_counts = new TYPE_INT_COUNTS_COLUMN[max_length + 1];

	int i,j;
	for (i= 0; i < PROFILEWIDTH; i++) 
		int_counts[0][i] = 0;

	int col = 1;

	while (col <= max_length) {
		if (input.eof() || input.peek() == EOF)
			break;

		input.read( (char*)&(int_counts[col]), PROFILEWIDTH * sizeof( int ) );
		col++;
	}
	col--;

	ImplProfile * p = new ImplProfile( regularizor, logoddor );

	p->setTrueLength(col);
	p->useFullLength();
	p->allocateCounts();
	CountColumn * counts = p->mCounts;

	for (i = 1; i <= col; i++) 
		for (j = 0; j < PROFILEWIDTH; j++) 
			counts[i][j] = (Count)(int_counts[i][j] / scale_factor);

	p->setPrepared(false);

	delete [] int_counts;

	return p;

}
//------------------------------------------------------------------------------------------
/** read counts of a profile from stream in binary format */
Alignandum * extractProfileBinaryCounts( std::istream & input, 
		const Position max_length,
		const Regularizor * regularizor,
		const LogOddor * logoddor ) {

	if (!regularizor)
		regularizor = getDefaultRegularizor();

	if (!logoddor)
		logoddor = getDefaultLogOddor();


	CountColumn * counts = new CountColumn[max_length + 1];
	for (int i= 0; i < PROFILEWIDTH; i++) 
		counts[0][i] = 0;

	int col = 1;

	while (col <= max_length) {
		if (input.eof() || input.peek() == EOF)
			break;

		input.read( (char*)&(counts[col]), PROFILEWIDTH * sizeof( Count) );
		col++;
	}

	ImplProfile * p = new ImplProfile( regularizor, logoddor );

	p->setTrueLength(col - 1);
	p->useFullLength();
	p->allocateCounts();

	memcpy( p->mCounts, counts, sizeof( CountColumn) * (col) );

	delete [] counts;

	p->setPrepared(false);

	return p;

}

//------------------------------------------------------------------------------------------
/** rescale counts from a profile by multiplying each entry by the scale_factor */
Alignandum * rescaleProfileCounts( Alignandum * dest,
		double scale_factor ) {

	// type cast to check, if we really have a profile
	ImplProfile * p_source = dynamic_cast<ImplProfile*>(dest);

	Position col, length;
	int i;
	length = p_source->getTrueLength();

	for ( col = 1; col <= length; col++) 
		for (i = 0; i < PROFILEWIDTH; i++) 
			p_source->mCounts[col][i] *= scale_factor;

	return dest;

}

//------------------------------------------------------------------------------------------
/** normalize counts from a profile so that all sum to total_weight per column*/
Alignandum * normalizeProfileCounts( Alignandum * dest,
		Count total_weight) {

	// type cast to check, if we really have a profile
	ImplProfile * p_source = dynamic_cast<ImplProfile*>(dest);

	Position col, length;
	int i;

	length = p_source->getTrueLength();

	for ( col = 1; col <= length; col++) {

		Count ntotal = 0;
		for (i = 0; i < PROFILEWIDTH; i++) 
			ntotal += p_source->mCounts[col][i];

		if (ntotal > 0) {
			double scale_factor = total_weight / ntotal;
			for (i = 0; i < PROFILEWIDTH; i++) 
				p_source->mCounts[col][i] *= scale_factor;
		}
	}

	return dest;

}

//------------------------------------------------------------------------------------------
/** substitutes columns in profile dest by columns in profile row using the mapping provided, where dest is in col and source is in row
 */
Alignandum * substituteProfileWithProfile( Alignandum * dest, const Alignandum * source, const Alignata * map_source2dest ) {

	// check, if we do have two profiles
	//!! to be implemented: some sensible warning messages
	const ImplProfile * p_source = dynamic_cast<const ImplProfile*>(source);
	const ImplProfile * p_dest = dynamic_cast<const ImplProfile*>(dest);

	AlignataConstIterator it(map_source2dest->begin());
	AlignataConstIterator it_end(map_source2dest->end());

	for (; it != it_end; ++it) {
		Position row = it->mRow;
		Position col = it->mCol;

		for (int i = 0; i < PROFILEWIDTH; i++) 
			p_dest->mCounts[col][i] = p_source->mCounts[row][i];
	}

	if (dest->isPrepared()) 
		dest->prepare();

	return dest;
}

//------------------------------------------------------------------------------------------
/** add counts of profile source to profile dest, using the mapping provided, where dest is in col and
    source is in row 
 */
Alignandum * addProfile2Profile( Alignandum * dest, const Alignandum * source, const Alignata * map_source2dest ) {

	// check, if we do have two profiles
	const ImplProfile * p_source = dynamic_cast<const ImplProfile*>(source);
	const ImplProfile * p_dest = dynamic_cast<const ImplProfile*>(dest);

	AlignataConstIterator it(map_source2dest->begin());
	AlignataConstIterator it_end(map_source2dest->end());

	for (; it != it_end; ++it) {
		Position row = it->mRow;
		Position col = it->mCol;

		for (int i = 0; i < PROFILEWIDTH; i++) 
			p_dest->mCounts[col][i] += p_source->mCounts[row][i];
	}

	if (dest->isPrepared()) {
		dest->release();	// first release, otherwise it won't calculate anew
		dest->prepare();
	}

	return dest;
}

//------------------------------------------------------------------------------------------
/** add sequence of source to profile dest, using the mapping provided, where dest is in col and
    source is in row 
 */
Alignandum * addSequence2Profile( Alignandum * dest, const Alignandum * source, const Alignata * map_source2dest ) {

	// check, if we do have two profiles
	const ImplProfile * p_dest = dynamic_cast<const ImplProfile*>(dest);

	AlignataConstIterator it(map_source2dest->begin());
	AlignataConstIterator it_end(map_source2dest->end());

	for (; it != it_end; ++it) {
		Position row = it->mRow;
		Position col = it->mCol;
		Residue r = source->asResidue(row);
		if (r < PROFILEWIDTH)
			p_dest->mCounts[col][r] ++;

	}

	if (dest->isPrepared()) {
		dest->release();	// first release, otherwise it won't calculate anew
		dest->prepare();
	}

	return dest;
}

//------------------------------------------------------------------------------------------
/** reset a profile to a new length. Clear old values.
 */
Alignandum * resetProfile( Alignandum * dest, Position new_length ) {

	// clear profile
	ImplProfile * profile = dynamic_cast<ImplProfile*>(dest);

	profile->release();
	profile->setTrueLength( new_length );
	profile->useFullLength();
	profile->allocateCounts();

	return dest;
}

//------------------------------------------------------------------------------------------
/** reset a profile to a new length. Clear old values.
 */
ProfileFrequencies * exportProfileFrequencies( Alignandum * dest ) {

	// clear profile
	ImplProfile * profile = dynamic_cast<ImplProfile*>(dest);

	Position from = dest->getFrom();
	Position to = dest->getTo();
	Position length = to - from + 1;

	bool was_prepared = false;

	if (!dest->isPrepared()) {
		dest->prepare();
		was_prepared = false;
	}

	ProfileFrequencies * result = new ProfileFrequencies(length + 1);
	unsigned int i = 1;

	(*result)[0].resize( PROFILEWIDTH, 0);

	for (Position col = from; col <= to; ++i, ++col) {
		(*result)[i].resize( PROFILEWIDTH, 0);
		for (unsigned int row = 0; row < PROFILEWIDTH; ++row) 
			(*result)[i][row] = profile->mFrequencies[col][row];
	}

	if (!was_prepared) 
		dest->release();

	return result;
}



//------------------------------------------------------------------------------------------
/** fill a profile with counts. The counts matrix has to be length + 1, the first row is
    not used
 */
/*
Alignandum * makeProfile( const CountsMatrix * src) {

  // check if counts are ok
  assert( src->getNumCols() == PROFILEWIDTH );

  // type cast to check, if we really have a profile
  ImplProfile * p_dest = dynamic_cast<ImplProfile*>(dest);

  // delete old counts and set new length
  p_dest->release();
  p_dest->setTrueLength( src->getNumRows() -1 );
  p_dest->allocateCounts();

  // retrieve pointer to member data
  AlignandumDataProfile & data   = (AlignandumDataProfile &)p_dest->getData();
  CountColumn  * profile_counts = data.mCountsPointer;

  // copy counts
  unsigned int row, col;

  for (row = 0; row < src->getNumRows(); row++)
    for (col = 0; col < PROFILEWIDTH; col++) 
      profile_counts[row][col] = src[row][col];

  return dest;
}

 */


} // namespace alignlib


