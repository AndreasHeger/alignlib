/*
  alignlib - a library for aligning protein sequences

  $Id: HelpersWeightor.cpp,v 1.2 2004/01/07 14:35:33 aheger Exp $

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
#include "Weightor.h"
#include "ImplWeightor.h"
#include "HelpersWeightor.h"

namespace alignlib {
  
  static const Weightor * DEFAULT_WEIGHTOR = makeNoWeightor();
  
  /** gets the default Weightor object */ 
  const Weightor * getDefaultWeightor() {
    return DEFAULT_WEIGHTOR;
  }

  /** sets the default Weightor object */
  const Weightor * setDefaultWeightor( const Weightor * weightor ) {
    const Weightor * t = DEFAULT_WEIGHTOR;
    DEFAULT_WEIGHTOR = weightor;
    return t;
  }

} // namespace alignlib
