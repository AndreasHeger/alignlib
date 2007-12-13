/*
  alignlib - a library for aligning protein sequences

  $Id: ImplProfile.cpp,v 1.3 2004/01/07 14:35:35 aheger Exp $

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
#include "ImplProfile.h"
#include "Translator.h"
#include "HelpersTranslator.h"
#include "LogOddor.h"
#include "Regularizor.h"
#include "MultipleAlignment.h"
#include "Alignata.h"
#include "AlignataIterator.h"
#include "HelpersAlignandum.h"
#include "HelpersLogOddor.h"
#include "HelpersRegularizor.h"

/** default objects */

#ifdef WITH_DMALLOC
#include <dmalloc.h>
#endif

using namespace std;

namespace alignlib {

//---------------------------------------> constructors and destructors <--------------------------------------
ImplProfile::ImplProfile( const Regularizor * regularizor, 
		const LogOddor * logoddor ) : 
			mRegularizor( regularizor ),
			mLogOddor( logoddor ),
			mCounts(NULL), mFrequencies(NULL), mProfile(NULL) 
			{

	if (regularizor == NULL) 
		mRegularizor = getDefaultRegularizor();

	if (logoddor == NULL)
		mLogOddor = getDefaultLogOddor();

			}

const AlignandumDataProfile & ImplProfile::getData() const 
{
	mData.mCountsPointer		= mCounts;
	mData.mFrequenciesPointer	= mFrequencies;
	mData.mProfilePointer		= mProfile;
	return mData;
}

//---------------------------------------------------------------------------------------------------------------
ImplProfile::ImplProfile(const ImplProfile & src ) : ImplAlignandum( src ), 
mRegularizor( src.mRegularizor ),
mLogOddor( src.mLogOddor),  
mCounts(NULL), 
mFrequencies(NULL), 
mProfile(NULL) 
{
	debug_func_cerr(5);

	// the first column is empty, so copy one more column
	Position copy_length = getFullLength(); 

	if (src.mCounts != NULL) {
		allocateCounts();
		memcpy( mCounts, src.mCounts, sizeof(CountColumn) * copy_length);
	}
	if (src.mFrequencies != NULL) {
		allocateFrequencies();
		memcpy( mFrequencies, src.mFrequencies, sizeof(FrequencyColumn) * copy_length);
	}

	if (src.mProfile != NULL) {
		allocateProfile();
		memcpy( mProfile, src.mProfile, sizeof(ProfileColumn) * copy_length);
	}

}

//---------------------------------------------------------------------------------------------------------------
ImplProfile::~ImplProfile() 
{
	debug_func_cerr(5);

	if (mCounts != NULL)
	{ delete [] mCounts; mCounts = NULL; }
	if (mFrequencies != NULL)
	{ delete [] mFrequencies; mFrequencies = NULL; }
	if (mProfile != NULL)
	{ delete [] mProfile; mProfile = NULL; }
}

//--------------------------------------------------------------------------------------
ImplProfile * ImplProfile::getClone() const {
	return new ImplProfile( *this );
}


//---------------------------------------------------------------------------------------------------------------
void ImplProfile::allocateCounts() 
{
	debug_func_cerr(5);

	if (mCounts != NULL)
		delete [] mCounts;
	mCounts = NULL;

	Position length = getFullLength();

	mCounts = new CountColumn[length];

	int i, j;
	for (i = 0; i < length; i++ )
		for (j = 0; j < PROFILEWIDTH; j++)
			mCounts[i][j] = 0;

}              

//---------------------------------------------------------------------------------------------------------------
void ImplProfile::allocateProfile() const 
{
	debug_func_cerr(5);

	if (mProfile != NULL)	
		delete [] mProfile;

	Position length = getFullLength();

	mProfile = new ProfileColumn[length];
	int i, j;
	for (i = 0; i < length; i++ )
		for (j = 0; j < PROFILEWIDTH; j++)
			mProfile[i][j] = 0;

}

//---------------------------------------------------------------------------------------------------------------
void ImplProfile::allocateFrequencies() const 
{
	debug_func_cerr(5);

	if (mFrequencies != NULL)
		delete [] mFrequencies;

	Position length = getFullLength();

	mFrequencies = new FrequencyColumn[length];

	int i, j;
	for (i = 0; i < length; i++ )
		for (j = 0; j < PROFILEWIDTH; j++)
			mFrequencies[i][j] = 0;
}

//---------------------------------------------------------------------------------------------------------------
Residue ImplProfile::getMaximumCounts( const Position column ) const 
{
	int i;

	Count max   = 0;
	Residue max_i   = getDefaultTranslator()->getGapCode();

	for (i = 0; i < PROFILEWIDTH; i++) {
		if (mCounts[column][i] > max) {
			max   = mCounts[column][i];
			max_i = i;
		}
	}  
	return max_i;
}

//---------------------------------------------------------------------------------------------------------------
Residue ImplProfile::getMaximumFrequencies( const Position column ) const {
	int i;

	Frequency max   = 0;
	Residue   max_i   = getDefaultTranslator()->getGapCode();

	for (i = 0; i < PROFILEWIDTH; i++) {
		if (mFrequencies[column][i] > max) {
			max   = mFrequencies[column][i];
			max_i = i;
		}
	}  
	return max_i;
}

//---------------------------------------------------------------------------------------------------------------
Residue ImplProfile::getMaximumProfile( const Position column ) const 
{
	int i;

	ProfileScore max   = 0;
	Residue   max_i   = getDefaultTranslator()->getGapCode();

	for (i = 0; i < PROFILEWIDTH; i++) {
		if (mProfile[column][i] > max) {
			max   = mProfile[column][i];
			max_i = i;
		}
	}  
	return max_i;
}

//---------------------------------------------------------------------------------------------------------------
void ImplProfile::mask( Position column) {
	int i;

	if (mCounts) 
		for (i = 0; i < PROFILEWIDTH; i++) 
			mCounts[column][i] = 0;

	if (mFrequencies)
		for (i = 0; i < PROFILEWIDTH; i++) 
			mFrequencies[column][i] = 0;

	if (mProfile) 
		for (i = 0; i < PROFILEWIDTH; i++) 
			mProfile[column][i] = MASK_VALUE;

}

//---------------------------------------------------------------------------------------------------------------
Residue ImplProfile::asResidue( const Position column ) const 
{

	if (mCounts) 
		return getMaximumCounts( column );

	return getDefaultTranslator()->getGapCode();
}

//--------------------------------------------------------------------------------------
void ImplProfile::prepare() const 
{
	debug_func_cerr(5);

	// do nothing, when a profile and frequencies already exist.
	if (!mFrequencies) {
		allocateFrequencies();
		mRegularizor->fillFrequencies( mFrequencies, mCounts, getFullLength() );
	}

	if (!mProfile) {
		allocateProfile();
		mLogOddor->fillProfile( mProfile, mFrequencies, getFullLength() );
	}
	setPrepared( true );
}



//--------------------------------------------------------------------------------------
void ImplProfile::release() const 
{
	if (mFrequencies != NULL)
	{
		delete [] mFrequencies;
		mFrequencies = NULL;
	}
	if (mProfile != NULL)
	{
		delete [] mProfile;
		mProfile = NULL;
	}
	setPrepared(false);
}

//--------------------------------------------------------------------------------------
void ImplProfile::shuffle( unsigned int num_iterations,
		Position window_size) {

	if (window_size == 0)
		window_size = getLength();

	Position first_from = getFrom();
	Count buffer[PROFILEWIDTH];

	for (unsigned x = 0; x < num_iterations; x++) { 

		Position i,j;
		Position to = getTo();

		while (to > first_from ) {
			Position from = to - window_size;

			if (from < 1) {
				from = 1;
				window_size = to;
			}

			for (i = to; i >= from; i--) {
				j = to - getRandomPosition(window_size) - 1;
				memcpy(buffer, mCounts[i], PROFILEWIDTH * sizeof(mCounts));
				memcpy(mCounts[i], mCounts[j], PROFILEWIDTH * sizeof(mCounts));
				memcpy(mCounts[j], buffer, PROFILEWIDTH * sizeof(mCounts));
			}

			to -= window_size;
		}
	}


	/*
  Position i,j;
  Position from = getFrom();
  Position from_1 = from - 1;
  Position to = getTo();
  Position length = getLength();

  Count buffer[PROFILEWIDTH];

  for (i = to; i >= from; i--) {
    j = GetRandomPosition(length) + from_1;
    memcpy(buffer, mCounts[i], PROFILEWIDTH * sizeof(mCounts));
    memcpy(mCounts[i], mCounts[j], PROFILEWIDTH * sizeof(mCounts));
    memcpy(mCounts[j], buffer, PROFILEWIDTH * sizeof(mCounts));
  }
	 */
	this->release();

}


//--------------------------------------> I/O <------------------------------------------------------------
void ImplProfile::write( std::ostream & output ) const {

	Position i;
	Position length = getLength();
	int j;

	output.setf( ios::fixed );

	if (mCounts) {
		output << "----------->counts<----------------------------------------" << endl;
		for (i = 0; i < length; i++) {
			output << setw(2) << i << "\t";
			for (j = 0; j < PROFILEWIDTH; j++) 
				output << setw(6) << setprecision(2) << mCounts[i][j];
			output << endl;

		}
	} else {
		output << "----------->no counts available<---------------------------" << endl;
	}

	if (mFrequencies) {
		output << "----------->frequencies<-----------------------------------" << endl;
		for (i = 0; i < length; i++) {
			output << setw(2) << i << "\t";
			for (j = 0; j < PROFILEWIDTH; j++) 
				output << setw(6) << setprecision(2) << mFrequencies[i][j];
			output << endl;
		}
	} else {
		output << "----------->no frequencies available<----------------------" << endl;
	}

	if (mProfile) {
		output << "----------->profile<---------------------------------------" << endl;
		for (i = 0; i < length; i++) {
			output << setw(2) << i << "\t";
			for (j = 0; j < PROFILEWIDTH; j++) 
				output << setw(6) << setprecision(2) << mProfile[i][j];
			output << endl;
		}
	} else {
		output << "----------->no profile available<--------------------------" << endl;
	}
}
//--------------------------------------------------------------------------------------
void ImplProfile::__save( std::ostream & output, MagicNumberType type ) const 
{
	if (type == MNNoType )
	{
		type = MNImplProfile;
		output.write( (char*)&type, sizeof(MagicNumberType ) );
	}		
	ImplAlignandum::__save( output, type );

	output.write( (char*)mCounts, sizeof( CountColumn) * getFullLength() );
	if (isPrepared() )
	{
		output.write( (char*)mFrequencies, sizeof( FrequencyColumn) * getFullLength() );
		output.write( (char*)mProfile, sizeof( ProfileColumn) * getFullLength() );	
	}		
}

//--------------------------------------------------------------------------------------
void ImplProfile::load( std::istream & input)  
{
	ImplAlignandum::load( input );

	allocateCounts();	
	input.read( (char*)mCounts, sizeof( CountColumn) * getFullLength() );

	if (input.fail() ) 
		throw AlignException( "incomplete profile in stream.");

	if (isPrepared() )
	{
		allocateFrequencies();
		input.read( (char*)mFrequencies, sizeof( FrequencyColumn) * getFullLength() );
		allocateProfile();
		input.read( (char*)mProfile, sizeof( ProfileColumn) * getFullLength() );	
	}		

}

//--------------------------------------> I/O <------------------------------------------------------------



} // namespace alignlib
