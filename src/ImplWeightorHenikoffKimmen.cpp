/*
  alignlib - a library for aligning protein sequences

  $Id: ImplWeightorHenikoffKimmen.cpp,v 1.2 2004/01/07 14:35:36 aheger Exp $

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
#include "alignlib_fwd.h"
#include "alignlib_interfaces.h"
#include "alignlib_fwd.h"
#include "AlignlibDebug.h"

#include "HelpersWeightor.h"
#include "Weightor.h"
#include "ImplWeightorHenikoffKimmen.h"
#include "MultipleAlignment.h"
#include "HelpersEncoder.h"
#include "Encoder.h"

using namespace std;

namespace alignlib 
{

/** factory functions */
HWeightor makeWeightorHenikoffKimmen() 
{ 
	return HWeightor(new ImplWeightorHenikoffKimmen());
}

//---------------------------------------------------------< constructors and destructors >--------------------------------------
ImplWeightorHenikoffKimmen::ImplWeightorHenikoffKimmen () : 
	ImplWeightorHenikoff() 
{
}

ImplWeightorHenikoffKimmen::~ImplWeightorHenikoffKimmen () 
{
}

ImplWeightorHenikoffKimmen::ImplWeightorHenikoffKimmen (const ImplWeightorHenikoffKimmen & src ) : 
	ImplWeightorHenikoff(src) 
{
}

//--------------------------------------------------------------------------------------------------------------------------------
HSequenceWeights ImplWeightorHenikoffKimmen::calculateWeights( 
		const HMultipleAlignment & src,
		const HEncoder & translator ) const 
{
	debug_func_cerr(5);

	HSequenceWeights weights( ImplWeightorHenikoff::calculateWeights( src, translator ) );

	int nsequences = src->getWidth();
	
	// rescale weights, so that they sum to the number of nsequences 
	rescaleWeights( weights, nsequences, nsequences);

	return weights;
}

} // namespace alignlib
