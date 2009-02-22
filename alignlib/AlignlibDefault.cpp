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
#include "HelpersToolkit.h"

using namespace std;

namespace alignlib
{

IMPLEMENT_STATIC_DEFAULT( HToolkit, makeToolkit(), getDefaultToolkit, setDefaultToolkit, default_toolkit );

// default objects without dependencies
IMPLEMENT_DEFAULT( HIterator2D, getDefaultIterator2D, setDefaultIterator2D, getIterator2D, setIterator2D );
IMPLEMENT_DEFAULT( HLogOddor, getDefaultLogOddor, setDefaultLogOddor,  getLogOddor, setLogOddor );
IMPLEMENT_DEFAULT( HRegularizor, getDefaultRegularizor, setDefaultRegularizor, getRegularizor, setRegularizor );
IMPLEMENT_DEFAULT( HWeightor, getDefaultWeightor, setDefaultWeightor, getWeightor, setWeightor );
IMPLEMENT_DEFAULT( HScorer, getDefaultScorer, setDefaultScorer,  getScorer, setScorer )
IMPLEMENT_DEFAULT( HDistor, getDefaultDistor, setDefaultDistor, getDistor, setDistor );
IMPLEMENT_DEFAULT( HEncoder, getDefaultEncoder, setDefaultEncoder, getEncoder, setEncoder );
IMPLEMENT_DEFAULT( HSubstitutionMatrix, getDefaultSubstitutionMatrix, setDefaultSubstitutionMatrix, getSubstitutionMatrix, setSubstitutionMatrix );
IMPLEMENT_DEFAULT( HTreetor, getDefaultTreetor, setDefaultTreetor, getTreetor, setTreetor );


IMPLEMENT_STATIC_DEFAULT( HPalette, makePalette(), getDefaultPalette, setDefaultPalette, default_palette );

} // namespace alignlib
