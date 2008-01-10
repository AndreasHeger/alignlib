/*
  alignlib - a library for aligning protein sequences

  $Id: HelpersRenderer.cpp,v 1.2 2004/01/07 14:35:32 aheger Exp $

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
#include "alignlib_fwd.h"
#include "alignlib_default.h"
#include "AlignlibDebug.h"
#include "ImplRenderer.h"
#include "HelpersRenderer.h"


using namespace std;

namespace alignlib 
{

const TYPE_PALETTE DEFAULT_PALETTE[16] = {
		/* default: Silver = */ "#C0C0C0", 
		/* Black	*/ "#000000", 
		/* Green	*/ "#008000",
		/* Lime	*/ "#00FF00",
		/* Gray	*/ "#808080", 
		/* Olive	*/ "#808000",
		/* White	*/ "#FFFFFF", 
		/* Yellow	*/ "#FFFF00",
		/* Maroon	*/ "#800000", 
		/* Navy	*/ "#000080",
		/* Red	*/ "#FF0000", 
		/* Blue	*/ "#0000FF",
		/* Purple	*/ "#800080", 
		/* Teal	*/ "#008080",
		/* Fuchsia	*/ "#FF00FF",
		/* Aqua	*/ "#00FFFF"
};

const TYPE_PALETTE * getDefaultPalette() { return DEFAULT_PALETTE; }

} // namespace alignlib
