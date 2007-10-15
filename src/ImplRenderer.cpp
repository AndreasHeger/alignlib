/*
  alignlib - a library for aligning protein sequences

  $Id: ImplRenderer.cpp,v 1.2 2004/01/07 14:35:36 aheger Exp $

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


#include "alignlib.h"
#include "AlignlibDebug.h"
#include "ImplRenderer.h"

#ifdef WITH_DMALLOC
#include <dmalloc.h>
#endif

using namespace std;

namespace alignlib {

//---------------------------------------------------------< constructors and destructors >--------------------------------------
ImplRenderer::ImplRenderer (const TYPE_PALETTE * palette) : mPalette ( palette ) {
}
		       
ImplRenderer::~ImplRenderer () {
}

ImplRenderer::ImplRenderer (const ImplRenderer & src ) : mPalette( src.mPalette ) {
}



} // namespace alignlib
