/*
  alignlib - a library for aligning protein sequences

  $Id: HelpersRegularizor.cpp,v 1.2 2004/01/07 14:35:32 aheger Exp $

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
#include "Regularizor.h"
#include "ImplRegularizor.h"
#include "HelpersRegularizor.h"

namespace alignlib {
  
  static Regularizor * DEFAULT_REGULARIZOR = makeRegularizorDirichlet();

  /** gets the default Regularizor object */ 
  const Regularizor * getDefaultRegularizor() {
    return DEFAULT_REGULARIZOR;
  }

  /** sets the default Regularizor object */
  void setDefaultRegularizor( Regularizor * regularizor ) 
  {
	  if (DEFAULT_REGULARIZOR != NULL)
		  DEFAULT_REGULARIZOR = regularizor;
	  DEFAULT_REGULARIZOR = regularizor;
  }

} // namespace alignlib
