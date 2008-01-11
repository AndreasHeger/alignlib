/*
  alignlib - a library for aligning protein sequences

  $Id: ImplAlignatorDummy.cpp,v 1.2 2004/01/07 14:35:34 aheger Exp $

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
#include "AlignlibDebug.h"

#include "ImplAlignatorDummy.h"

using namespace std;

namespace alignlib 
{

HAlignator makeAlignatorDummy( const HAlignment & ali)
{
	return HAlignator( new ImplAlignatorDummy( ali ) );
}

//---------------------------------------------------------< constructors and destructors >--------------------------------------
ImplAlignatorDummy::ImplAlignatorDummy ( const HAlignment & ali) : 
	ImplAlignator(), mAlignment ( ali )
	{
	}

ImplAlignatorDummy::~ImplAlignatorDummy ()
{
}

ImplAlignatorDummy::ImplAlignatorDummy (const ImplAlignatorDummy & src ) :
	ImplAlignator(src), 
	mAlignment(src.mAlignment)
	{
	}

//----------------------------------------------------------------------------------------------------------
HAlignator ImplAlignatorDummy::getClone() const 
{
	return HAlignator( new ImplAlignatorDummy( *this ) );
}

//----------------------------------------------------------------------------------------------------------
HAlignment & ImplAlignatorDummy::align( HAlignment & result,
		const HAlignandum & row, 
		const HAlignandum & col ) 
{

	startUp(result, row, col );

	// copy alignment only in region given and move to (1,1)
	copyAlignment( result, mAlignment,
			mIterator->row_front(), mIterator->row_back(),
			mIterator->col_front(), mIterator->col_back() );

	cleanUp( result, row, col );

	return result;
}


} // namespace alignlib
