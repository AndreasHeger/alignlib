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

#include "alignlib_fwd.h"
#include "alignlib_interfaces.h"
#include "AlignlibDebug.h"
#include "AlignException.h"
#include "ImplProfile.h"
#include "Encoder.h"
#include "HelpersEncoder.h"
#include "LogOddor.h"
#include "Weightor.h"
#include "Regularizor.h"
#include "MultipleAlignment.h"
#include "Alignment.h"
#include "AlignmentIterator.h"
#include "HelpersAlignandum.h"
#include "HelpersLogOddor.h"
#include "HelpersRegularizor.h"
#include "HelpersMultipleAlignment.h"
#include "HelpersWeightor.h"

/** default objects */

using namespace std;

namespace alignlib 
{

//---------------------------------< implementation of factory functions >--------------

//------------------------------------------------------------------------------------------
/** create empty profile 
 * */
HAlignandum makeProfile( 
		const HEncoder & translator,
		const HWeightor & weightor,
		const HRegularizor & regularizor,
		const HLogOddor & logoddor ) 
{
	return HAlignandum( new ImplProfile( translator, weightor, regularizor, logoddor ) );
}

HAlignandum makeProfile()
{
	return makeProfile( 
			getDefaultEncoder(),
			getDefaultWeightor(),
			getDefaultRegularizor(),
			getDefaultLogOddor() );
}

//------------------------------------------------------------------------------------------
/** create empty profile with given length */
HAlignandum makeProfile( const Position & length,
		const HEncoder & translator,
		const HWeightor & weightor,		
		const HRegularizor & regularizor,
		const HLogOddor & logoddor ) 
		{

	return HAlignandum( new ImplProfile( length, 
			translator, weightor, regularizor, logoddor ) );
}

HAlignandum makeProfile( const Position & length )
{
	return makeProfile( length, 
			getDefaultEncoder(),
			getDefaultWeightor(),
			getDefaultRegularizor(),
			getDefaultLogOddor() );
}


//------------------------------------------------------------------------------------------
/** create profile from a string of sequences */
HAlignandum makeProfile( const std::string & src, 
		int nsequences,
		const HEncoder & translator,
		const HWeightor & weightor, 
		const HRegularizor & regularizor,
		const HLogOddor & logoddor ) 
		{

	HMultipleAlignment m( makeMultipleAlignment() );
	fillMultipleAlignment( 
			m,
			src, 
			nsequences );

	return HAlignandum (new ImplProfile(
			m, 
			translator,
			weightor,
			regularizor, 
			logoddor ));
}

HAlignandum makeProfile( const std::string & src, int nsequences)
{
	return makeProfile( src, nsequences,
			getDefaultEncoder(),
			getDefaultWeightor(),
			getDefaultRegularizor(),
			getDefaultLogOddor() );
}


//------------------------------------------------------------------------------------------
/** create a default profile from a multiple alignment */
HAlignandum makeProfile( 
		const HMultipleAlignment & mali, 
		const HEncoder & translator,
		const HWeightor & weightor, 
		const HRegularizor & regularizor,
		const HLogOddor & logoddor ) 
		{

	return HAlignandum( new ImplProfile( mali, 
			translator, 
			weightor,
			regularizor, 
			logoddor ));
}

HAlignandum makeProfile( const HMultipleAlignment & mali )
{
	return makeProfile( mali,
			getDefaultEncoder(), 
			getDefaultWeightor(),
			getDefaultRegularizor(),
			getDefaultLogOddor());
}

//---------------------------------------> constructors and destructors <--------------------------------------
// The constructor is potentially empty, so that this object can be read from file.
ImplProfile::ImplProfile( 
		const HEncoder & translator,
		const HWeightor & weightor,
		const HRegularizor & regularizor, 
		const HLogOddor & logoddor ) :
			ImplAlignandum( translator ),
			mWeightor( weightor ),
			mRegularizor( regularizor ),
			mLogOddor( logoddor ),
			mCountMatrix(NULL), mFrequencyMatrix(NULL), mScoreMatrix(NULL),
			mProfileWidth(0)
{
	debug_func_cerr(5);
	allocateCounts();
}

ImplProfile::ImplProfile(
		const Position & length,
		const HEncoder & translator,
		const HWeightor & weightor,
		const HRegularizor & regularizor, 
		const HLogOddor & logoddor ) :
			ImplAlignandum( translator ),
			mWeightor( weightor ),
			mRegularizor( regularizor ),
			mLogOddor( logoddor ),
			mCountMatrix(NULL), mFrequencyMatrix(NULL), mScoreMatrix(NULL),
			mProfileWidth(0)
{
	debug_func_cerr(5);
	resize( length );
}


ImplProfile::ImplProfile( 
		const HMultipleAlignment & src,
		const HEncoder & translator,
		const HWeightor & weightor,
		const HRegularizor & regularizor, 
		const HLogOddor & logoddor ) :
			ImplAlignandum( translator ),
			mWeightor( weightor ),
			mRegularizor( regularizor ),
			mLogOddor( logoddor ),
			mCountMatrix(NULL), mFrequencyMatrix(NULL), mScoreMatrix(NULL),
			mProfileWidth(0)
{
	debug_func_cerr(5);
	resize( src->getLength() );
	fillCounts( src );
}

		
		
//---------------------------------------------------------------------------------------------------------------
ImplProfile::ImplProfile(const ImplProfile & src ) : ImplAlignandum( src ), 
	mRegularizor( src.mRegularizor ),
	mLogOddor( src.mLogOddor),  
	mCountMatrix(NULL), 
	mFrequencyMatrix(NULL), 
	mScoreMatrix(NULL),
	mProfileWidth(src.mProfileWidth)
{
	debug_func_cerr(5);

	if (src.mCountMatrix != NULL) 
		mCountMatrix = new CountMatrix( *src.mCountMatrix );		

	if (src.mFrequencyMatrix != NULL) 
		mFrequencyMatrix = new FrequencyMatrix( *src.mFrequencyMatrix );
	
	if (src.mScoreMatrix != NULL) 
		mScoreMatrix = new ScoreMatrix( *src.mScoreMatrix );
}

//---------------------------------------------------------------------------------------------------------------
ImplProfile::~ImplProfile() 
{
	debug_func_cerr(5);

	if (mCountMatrix != NULL)
		{ delete mCountMatrix; mCountMatrix = NULL; }
	if (mFrequencyMatrix != NULL)
		{ delete mFrequencyMatrix; mFrequencyMatrix = NULL; }
	if (mScoreMatrix != NULL)
		{ delete mScoreMatrix; mScoreMatrix = NULL; }
}

//--------------------------------------------------------------------------------------
HAlignandum ImplProfile::getClone() const 
{
	return HAlignandum( new ImplProfile( *this ) );
}

//--------------------------------------------------------------------------------------
CountMatrix * ImplProfile::getCountMatrix() const
{
	return mCountMatrix;
}

//--------------------------------------------------------------------------------------
FrequencyMatrix * ImplProfile::getFrequencyMatrix() const
{
	return mFrequencyMatrix;
}

//--------------------------------------------------------------------------------------
ScoreMatrix * ImplProfile::getScoreMatrix() const
{
	return mScoreMatrix;
}

//--------------------------------------------------------------------------------------
void ImplProfile::resize( Position length )
{
	ImplAlignandum::resize( length );
	allocateCounts();
}

//---------------------------------------------------------------------------------------------------------------
void ImplProfile::fillCounts( const HMultipleAlignment &  src )
{
	debug_func_cerr(5);

	resize( src->getLength() );
	mWeightor->fillCounts( *mCountMatrix, src, mEncoder );

	setPrepared( false );  
}

//---------------------------------------------------------------------------------------------------------------
template< class T>
Matrix<T> * ImplProfile::allocateSegment( Matrix<T> * data ) const
{
	debug_func_cerr(5);

	if (data != NULL) delete data;
	data = new Matrix<T>( getFullLength(), mProfileWidth, 0 );
	
	return data;
}

//---------------------------------------------------------------------------------------------------------------
void ImplProfile::allocateCounts() const
{
	debug_func_cerr(5);
	
	debug_cerr( 5, "allocating counts for a profile of size " << getFullLength() << " x " << (int)mProfileWidth );
	
	mProfileWidth = mEncoder->getAlphabetSize();
	mCountMatrix = allocateSegment<Count>( mCountMatrix );
}              

//---------------------------------------------------------------------------------------------------------------
void ImplProfile::allocateScores() const 
{
	debug_func_cerr(5);
	mScoreMatrix = allocateSegment<Score>( mScoreMatrix );
}

//---------------------------------------------------------------------------------------------------------------
void ImplProfile::allocateFrequencies() const 
{
	debug_func_cerr(5);
	mFrequencyMatrix = allocateSegment<Frequency>( mFrequencyMatrix );
}

//---------------------------------------------------------------------------------------------------------------
template<class T>
Residue ImplProfile::getMaximumPerColumn( const Matrix<T> * data, 
										  const Position & column ) const
{
	assert( data != NULL );
	
	T max = std::numeric_limits<T>::min();
	Residue max_i = 0;

	T * col = data->getRow( column );
	for (Residue i = 0; i < mProfileWidth; ++i) 
	{
		if (col[i] > max) 
		{
			max = col[i];
			max_i = i;
		}
	}  
	
	// if no counts, return gap code
	if (max == 0)
		return mEncoder->getGapCode();
	
	return max_i;
}

//---------------------------------------------------------------------------------------------------------------
Residue ImplProfile::getMaximumCount( const Position column ) const 
{
	debug_func_cerr(6);
	return getMaximumPerColumn<Count>( mCountMatrix, column );
}

//---------------------------------------------------------------------------------------------------------------
Residue ImplProfile::getMaximumFrequency( const Position column ) const 
{
	debug_func_cerr(5);
	return getMaximumPerColumn<Frequency>( mFrequencyMatrix, column );
}

//---------------------------------------------------------------------------------------------------------------
Residue ImplProfile::getMaximumScore( const Position column ) const 
{
	debug_func_cerr(5);	
	return getMaximumPerColumn<Score>( mScoreMatrix, column );
}

//---------------------------------------------------------------------------------------------------------------
template<class T>
void ImplProfile::setColumnToValue( Matrix<T> * data, 
		const Position & column,
		const T & value )
{
	if (data == NULL) return;
		
	T * col = data->getRow( column );
	for (int i = 0; i < mProfileWidth; i++) 
		col[i] = 0;
}

//---------------------------------------------------------------------------------------------------------------
void ImplProfile::mask( const Position & column) 
{
	ImplAlignandum::mask( column );
	setColumnToValue<Count>( mCountMatrix, column, 0);
	setColumnToValue<Frequency>( mFrequencyMatrix, column, 0);
	setColumnToValue<Score>( mScoreMatrix, column, 0 );
}

//---------------------------------------------------------------------------------------------------------------
Residue ImplProfile::asResidue( const Position column ) const 
{
	if (isMasked(column))
		return mEncoder->getMaskCode();
	return getMaximumCount( column );
}

//--------------------------------------------------------------------------------------
void ImplProfile::prepare() const 
{
	debug_func_cerr(5);

	// do nothing, when a profile and frequencies already exist.
	if (mFrequencyMatrix == NULL) 
	{
		allocateFrequencies();
		mRegularizor->fillFrequencies( *mFrequencyMatrix, *mCountMatrix, mEncoder ); 
	}

	if (!mScoreMatrix) 
	{
		allocateScores();
		mLogOddor->fillProfile( *mScoreMatrix, *mFrequencyMatrix, mEncoder );
	}
	setPrepared( true );
}

//--------------------------------------------------------------------------------------
void ImplProfile::release() const 
{
	debug_func_cerr(5);
	
	if (mFrequencyMatrix != NULL)
	{
		delete mFrequencyMatrix;
		mFrequencyMatrix = NULL;
	}
	if (mScoreMatrix != NULL)
	{
		delete mScoreMatrix;
		mScoreMatrix = NULL;
	}
	setPrepared(false);
}

//--------------------------------------------------------------------------------------
void ImplProfile::swap( const Position & x, const Position & y )
{
	mCountMatrix->swapRows( x, y );
	if (mFrequencyMatrix != NULL)
		mFrequencyMatrix->swapRows( x, y );
	if (mCountMatrix != NULL)
		mCountMatrix->swapRows( x, y );
}

//--------------------------------------------------------------------------------------
template<class T>
void ImplProfile::writeSegment( std::ostream & output, const Matrix<T> * data ) const
{
	if (data == NULL) return;
	
	output << setw(5) << "#" << "  " << " ";
	for (Residue j = 0; j < mProfileWidth; j++) 
		output << setw(6) << mEncoder->decode( j );
	output << std::endl;
	for (int i = 0; i < getLength(); i++) 
	{
		output << setw(5) << i << " " << asChar(i) << " ";
		const T * column = data->getRow( i );		
		for (Residue j = 0; j < mProfileWidth; j++) 
			output << setw(6) << setprecision(2) << column[j];
		output << endl;
	}
}

//--------------------------------------> I/O <------------------------------------------------------------
void ImplProfile::write( std::ostream & output ) const 
{

	output.setf( ios::fixed );

	if (mCountMatrix) 
	{
		output << "----------->counts<----------------------------------------" << endl;
		writeSegment<Count>( output, mCountMatrix );
	}
	else 
	{
		output << "----------->no counts available<---------------------------" << endl;
	}

	if (mFrequencyMatrix) 
	{
		output << "----------->frequencies<-----------------------------------" << endl;
		writeSegment<Frequency>( output, mFrequencyMatrix );
	}
	else 
	{
		output << "----------->no frequencies available<----------------------" << endl;
	}

	if (mScoreMatrix) 
	{
		output << "----------->profile<---------------------------------------" << endl;
		writeSegment<Score>( output, mScoreMatrix );		
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

	size_t size = getFullLength() * mProfileWidth;
	
	output.write( (char*)mCountMatrix->getData(), sizeof(Count) * size);
	if (isPrepared() )
	{
		output.write( (char*)mFrequencyMatrix->getData(), sizeof(Frequency) * size);
		output.write( (char*)mScoreMatrix->getData(), sizeof(Score) * size );	
	}		
}

//--------------------------------------------------------------------------------------
void ImplProfile::load( std::istream & input)  
{
	ImplAlignandum::load( input );

	input.read( (char*)&mProfileWidth, sizeof( Residue ) );
	
	allocateCounts();
	
	size_t size = getFullLength() * mProfileWidth;
	input.read( (char*)mCountMatrix->getData(), 
			sizeof( Count) * size );

	if (input.fail() ) 
		throw AlignException( "incomplete profile in stream.");

	if (isPrepared() )
	{
		allocateFrequencies();
		input.read( (char*)mFrequencyMatrix->getData(), 
				sizeof( Frequency) * size );
		allocateScores();
		input.read( (char*)mScoreMatrix->getData(), 
				sizeof(Score) * size );
	}
}

//--------------------------------------> I/O <------------------------------------------------------------

} // namespace alignlib
