/*
  alignlib - a library for aligning protein sequences

  $Id: ImplLogOddorBackground.cpp,v 1.2 2004/01/07 14:35:35 aheger Exp $

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
#include <math.h>
#include "alignlib.h"
#include "alignlib_fwd.h"
#include "AlignException.h"
#include "AlignlibDebug.h"
#include "ImplLogOddorBackground.h"
#include "HelpersLogOddor.h"

using namespace std;

namespace alignlib 
{

//---------------------------------------------------------< factory functions >--------------------------------------
HLogOddor makeLogOddorBackground( const FrequencyVector & frequencies,
		const Score & scale, const Score & mask_value )
{
	return HLogOddor(new ImplLogOddorBackground( frequencies, scale, mask_value ));
}

//---------------------------------------------------------< constructors and destructors >--------------------------------------
ImplLogOddorBackground::ImplLogOddorBackground ( const FrequencyVector & frequencies,
		const Score & scale_factor, const Score & mask_value ) :
	ImplLogOddor( scale_factor, mask_value ),
	mBackgroundFrequencies( frequencies )
	{
	}

ImplLogOddorBackground::~ImplLogOddorBackground () 
{
}

ImplLogOddorBackground::ImplLogOddorBackground (const ImplLogOddorBackground & src ) :
	ImplLogOddor( src ), mBackgroundFrequencies( src.mBackgroundFrequencies )
	{
	}

//--------------------------------------------------------------------------------------------------------------------------------
void ImplLogOddorBackground::fillProfile( ScoreMatrix * profile ,
		const FrequencyMatrix * frequencies ) const 
		{
	debug_func_cerr(5);
	
	// simply take the frequencies and divide by background-frequencies and take log. 
	// For frequencies of 0, MASK_VALUE is used.
	Position length = frequencies->getNumRows();
	Residue width  = frequencies->getNumCols();

	if (mBackgroundFrequencies.size() != width )
		throw AlignException("ImplLogOddorBackground: the size of alphabet does not correspond to number of supplied background frequencies.");

	for (Position column = 0; column < length; column++) 
	{
		const Frequency * fcolumn = frequencies->getRow(column);
		Score * pcolumn = profile->getRow(column);
		for (Residue i = 0; i < width; ++i)
		{
			Frequency f = 0;	
			if ((f = fcolumn[i]) > 0)
				pcolumn[i] = log(f / mBackgroundFrequencies[i]) / mScaleFactor;
			else
				pcolumn[i] = mMaskValue;
		}
	}
		}

} // namespace alignlib








