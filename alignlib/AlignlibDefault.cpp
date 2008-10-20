/*
  alignlib - a library for aligning protein sequences

  $Id: Alignment.cpp,v 1.2 2004/01/07 14:35:31 aheger Exp $

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
#include "alignlib_fwd.h"
#include "alignlib_interfaces.h"
#include "alignlib_default.h"

#include "HelpersEncoder.h"
#include "HelpersDistor.h"
#include "HelpersRegularizor.h"
#include "HelpersPalette.h"
#include "HelpersIterator2D.h"
#include "HelpersScorer.h"
#include "HelpersLogOddor.h"
#include "HelpersDistor.h"
#include "HelpersWeightor.h"
#include "HelpersSubstitutionMatrix.h"
#include "HelpersTreetor.h"

using namespace std;

namespace alignlib 
{

// default objects without dependencies
IMPLEMENT_DEFAULT( HIterator2D, makeIterator2DFull(), getDefaultIterator2D, setDefaultIterator2D, default_iterator2d );
IMPLEMENT_DEFAULT( HLogOddor, makeLogOddor(), getDefaultLogOddor, setDefaultLogOddor, default_logoddor );
IMPLEMENT_DEFAULT( HRegularizor, makeRegularizor(), getDefaultRegularizor, setDefaultRegularizor, default_regularizor );
IMPLEMENT_DEFAULT( HWeightor, makeWeightor(), getDefaultWeightor, setDefaultWeightor, default_weightor );
IMPLEMENT_DEFAULT( HPalette, makePalette(), getDefaultPalette, setDefaultPalette, default_palette );
IMPLEMENT_DEFAULT( HScorer, makeScorer(), getDefaultScorer, setDefaultScorer, default_scorer )
IMPLEMENT_DEFAULT( HDistor, makeDistorClustal(), getDefaultDistor, setDefaultDistor, default_distor );

// this depends on the default alphabets already being contstructed
// If problems, switch to makeEncoder idiom like for other objects.
IMPLEMENT_DEFAULT( HEncoder, getEncoder( Protein20 ), 
		getDefaultEncoder, setDefaultEncoder, default_translator );

// default objects with dependencies
IMPLEMENT_DEFAULT( HSubstitutionMatrix, 
		makeSubstitutionMatrixBlosum62( getDefaultEncoder() ), 
		getDefaultSubstitutionMatrix,
		setDefaultSubstitutionMatrix,
		default_matrix )

IMPLEMENT_DEFAULT( HTreetor, 
		makeTreetorDistanceLinkage( getDefaultDistor() ), 
		getDefaultTreetor, 
		setDefaultTreetor,
		default_treetor );

} // namespace alignlib