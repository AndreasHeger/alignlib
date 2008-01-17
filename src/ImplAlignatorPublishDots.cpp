/*
  alignlib - a library for aligning protein sequences

  $Id: ImplAlignatorPublishDots.cpp,v 1.2 2004/01/07 14:35:34 aheger Exp $

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
#include "alignlib_fwd.h"
#include "alignlib_interfaces.h"
#include "AlignlibDebug.h"
#include "Alignandum.h"

#include "Alignment.h"
#include "HelpersAlignment.h"

#include "SubstitutionMatrix.h"
#include "HelpersSubstitutionMatrix.h"

#include "ImplAlignatorPublishAlignment.h"

using namespace std;

namespace alignlib 
{

HAlignator makeAlignatorPublishAlignment( HAlignment ali) {
  return HAlignator(new ImplAlignatorPublishAlignment( ali ) );
}

//---------------------------------------------------------< constructors and destructors >--------------------------------------
ImplAlignatorPublishAlignment::ImplAlignatorPublishAlignment ( HAlignment ali) : 
  ImplAlignator( getDefaultSubstitutionMatrix(), 0, 0 ), 
  mAlignment ( ali ) 
  {
}
		       
ImplAlignatorPublishAlignment::~ImplAlignatorPublishAlignment () 
{
}

ImplAlignatorPublishAlignment::ImplAlignatorPublishAlignment (const ImplAlignatorPublishAlignment & src ) : ImplAlignator(src) 
{
}

//--------------------------------------------------------------------------------------------------------
ImplAlignatorPublishDots * ImplAlignatorPublishDots::getClone() const 
{
 return new ImplAlignatorPublishDots( *this );
}


void ImplAlignatorPublishAlignment::align( const HAlignandum row, const HAlignandum col, HAlignment result) 
{
  startUp(row, col, result );
  
  cleanUp( row, col, results );
}


} // namespace alignlib
