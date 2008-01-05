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
#include "Alignment.h"
#include "AlignmentIterator.h"
#include "HelpersAlignandum.h"
#include "HelpersLogOddor.h"
#include "HelpersRegularizor.h"

/** default objects */

using namespace std;

namespace alignlib 
{

//---------------------------------------> constructors and destructors <--------------------------------------
ImplProfile::ImplProfile( const Translator * translator, 
		const Regularizor * regularizor, 
		const LogOddor * logoddor ) :
			ImplAlignandum( translator ),
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

	if (src.mCounts != NULL) 
	{
		allocateCounts();
		memcpy( mCounts, src.mCounts, sizeof(Count) * mProfileWidth * copy_length);
	}
	if (src.mFrequencies != NULL) 
	{
		allocateFrequencies();
		memcpy( mFrequencies, src.mFrequencies, sizeof(Frequency) * mProfileWidth * copy_length);
	}

	if (src.mProfile != NULL) 
	{
		allocateProfile();
		memcpy( mProfile, src.mProfile, sizeof(Score) * mProfileWidth * copy_length);
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
ImplProfile * ImplProfile::getClone() const 
{
	return new ImplProfile( *this );
}

//---------------------------------------------------------------------------------------------------------------
template< class T>
void ImplProfile::allocateSegment( T * data ) const
{
	debug_func_cerr(5);

	if (data != NULL)
		delete [] data;
	data = NULL;

	Position length = getFullLength();

	data = new T[length * mProfileWidth];

	int i, j;
	for (i = 0; i < length * mProfileWidth; i++ )
		data[i] = 0;
}

//---------------------------------------------------------------------------------------------------------------
void ImplProfile::allocateCounts() const
{
	debug_func_cerr(5);
	allocateSegment<Count>( mCounts );
}              

//---------------------------------------------------------------------------------------------------------------
void ImplProfile::allocateProfile() const 
{
	debug_func_cerr(5);
	allocateSegment<Score>( mProfile );
}

//---------------------------------------------------------------------------------------------------------------
void ImplProfile::allocateFrequencies() const 
{
	debug_func_cerr(5);
	allocateSegment<Score>( mProfile );
}

//---------------------------------------------------------------------------------------------------------------
template<class T>
Residue ImplProfile::getMaximumPerColumn( const T * data, 
										  const Position & column ) const
{
	Count max = 0;
	Residue max_i = getDefaultTranslator()->getGapCode();

	Count * col = &mCounts[column * mProfileWidth];
	for (int i = 0; i < mProfileWidth; i++) 
	{
		if (col[i] > max) 
		{
			max   = col[i];
			max_i = i;
		}
	}  
	return max_i;
}

//---------------------------------------------------------------------------------------------------------------
Residue ImplProfile::getMaximumCounts( const Position column ) const 
{
	return getMaximumPerColumn<Count>( mCounts, column );
}

//---------------------------------------------------------------------------------------------------------------
Residue ImplProfile::getMaximumFrequencies( const Position column ) const 
{
	return getMaximumPerColumn<Frequency>( mFrequencies, column );
}

//---------------------------------------------------------------------------------------------------------------
Residue ImplProfile::getMaximumProfile( const Position column ) const 
{
	return getMaximumPerColumn<Score>( mProfile, column );
}

//---------------------------------------------------------------------------------------------------------------
template<class T>
void ImplProfile::maskColumn( T * data, const Position column )
{
	if (data == NULL) return;
		
	T * col = &data[column * mProfileWidth ];
	for (int i = 0; i < mProfileWidth; i++) 
		col[i] = 0;
}

//---------------------------------------------------------------------------------------------------------------
void ImplProfile::mask( const Position column) 
{
	maskColumn<Count>( mCounts, column);
	maskColumn<Frequency>( mFrequencies, column);
	maskColumn<Score>( mProfile, column);
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
	if (!mFrequencies) 
	{
		allocateFrequencies();
		mRegularizor->fillFrequencies( mFrequencies, 
				mCounts, 
				getFullLength(), 
				mProfileWidth );
	}

	if (!mProfile) 
	{
		allocateProfile();
		mLogOddor->fillProfile( mProfile, 
				mFrequencies, 
				getFullLength(),
				mProfileWidth); 
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
		Position window_size) 
{

	if (window_size == 0)
		window_size = getLength();

	Position first_from = getFrom();
	Count buffer[mProfileWidth];

	for (unsigned x = 0; x < num_iterations; x++) 
	{ 

		Position i,j;
		Position to = getTo();

		while (to > first_from ) 
		{
			Position from = to - window_size;

			if (from < 1) 
			{
				from = 1;
				window_size = to;
			}

			for (i = to; i >= from; i--) 
			{
				j = to - getRandomPosition(window_size) - 1;
				memcpy(buffer, 
						&mCounts[i * mProfileWidth], 
						mProfileWidth * sizeof(mCounts));
				memcpy(&mCounts[i*mProfileWidth], 
						&mCounts[j*mProfileWidth], 
						mProfileWidth * sizeof(mCounts));
				memcpy(&mCounts[j*mProfileWidth], 
						buffer, 
						mProfileWidth * sizeof(mCounts));
			}

			to -= window_size;
		}
	}
	this->release();
}

//--------------------------------------------------------------------------------------
template<class T>
void ImplProfile::writeSegment( std::ostream & output, const T * data ) const
{
	if (data == NULL) return;
	
	for (int i = 0; i < getLength(); i++) 
	{
		output << setw(2) << i << "\t";
		const T * column = &data[ i * mProfileWidth ];		
		for (Residue j = 0; j < mProfileWidth; j++) 
			output << setw(6) << setprecision(2) << column[j];
		output << endl;
	}
}

//--------------------------------------> I/O <------------------------------------------------------------
void ImplProfile::write( std::ostream & output ) const 
{

	output.setf( ios::fixed );

	if (mCounts) 
	{
		output << "----------->counts<----------------------------------------" << endl;
		writeSegment<Count>( output, mCounts );
	}
	else 
	{
		output << "----------->no counts available<---------------------------" << endl;
	}

	if (mFrequencies) 
	{
		output << "----------->frequencies<-----------------------------------" << endl;
		writeSegment<Frequency>( output, mFrequencies );
	}
	else 
	{
		output << "----------->no frequencies available<----------------------" << endl;
	}

	if (mProfile) 
	{
		output << "----------->profile<---------------------------------------" << endl;
		writeSegment<Score>( output, mProfile );		
	} 
	else 
	{
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

	output.write( (char*)&mProfileWidth, sizeof(Residue) );
	
	output.write( (char*)mCounts, sizeof(Count) * getFullLength() * mProfileWidth );
	if (isPrepared() )
	{
		output.write( (char*)mFrequencies, sizeof(Frequency) * getFullLength() * mProfileWidth );
		output.write( (char*)mProfile, sizeof(Score) * getFullLength() * mProfileWidth );	
	}		
}

//--------------------------------------------------------------------------------------
void ImplProfile::load( std::istream & input)  
{
	ImplAlignandum::load( input );

	input.read( (char*)&mProfileWidth, sizeof( Residue ) );
	
	allocateCounts();	
	input.read( (char*)mCounts, sizeof( Count) * getFullLength() * mProfileWidth);

	if (input.fail() ) 
		throw AlignException( "incomplete profile in stream.");

	if (isPrepared() )
	{
		allocateFrequencies();
		input.read( (char*)mFrequencies, 
				sizeof( Frequency) * getFullLength() * mProfileWidth );
		allocateProfile();
		input.read( (char*)mProfile, 
				sizeof(Score) * getFullLength() * mProfileWidth );	
	}		
}

//--------------------------------------> I/O <------------------------------------------------------------

} // namespace alignlib
