/*
  alignlib - a library for aligning protein sequences

  $Id: ImplLogOddor.cpp,v 1.2 2004/01/07 14:35:35 aheger Exp $

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
#include "AlignlibDebug.h"
#include "ImplLogOddor.h"
#include "HelpersLogOddor.h"

#ifdef WITH_DMALLOC
#include <dmalloc.h>
#endif

using namespace std;

namespace alignlib {

//---------------------------------------------------------< constructors and destructors >--------------------------------------
ImplLogOddor::ImplLogOddor ( const Frequency * f, Score scale_factor) : 
  mBackgroundFrequencies( f ), mScaleFactor( scale_factor ){
}
		       
ImplLogOddor::~ImplLogOddor () {
}

ImplLogOddor::ImplLogOddor (const ImplLogOddor & src ) : mScaleFactor (src.mScaleFactor ) {
}


//--------------------------------------------------------------------------------------------------------------------------------
void ImplLogOddor::fillProfile( ProfileColumn * profile ,
				const FrequencyColumn * frequencies, 
				const Position length ) const 
{
  debug_func_cerr(5);


  // simply take the frequencies and divide by background-frequencies and take log. 
  // For frequencies of 0, MASK_VALUE is used.
 
  Position column;
  Frequency f;

  for (column = 0; column < length; column++) 
      for (int i = 0; i < PROFILEWIDTH; i++) 
	if ((f = frequencies[column][i]) > 0)
	  profile[column][i] = log(f / mBackgroundFrequencies[i]) / mScaleFactor;
	else
	  profile[column][i] = MASK_VALUE;
}

} // namespace alignlib








