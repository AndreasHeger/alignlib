/*
  alignlib - a library for aligning protein sequences

  $Id: ImplNoWeightor.cpp,v 1.2 2004/01/07 14:35:35 aheger Exp $

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
#include "alignlib.h"
#include "AlignlibDebug.h"

#include "Weightor.h"
#include "ImplNoWeightor.h"
#include "MultipleAlignment.h"

#include "HelpersTranslator.h"

#ifdef WITH_DMALLOC
#include <dmalloc.h>
#endif

using namespace std;

namespace alignlib {


    /** factory functions */
    Weightor * makeNoWeightor() { return new ImplNoWeightor( getDefaultTranslator() );}

  
//---------------------------------------------------------< constructors and destructors >--------------------------------------
ImplNoWeightor::ImplNoWeightor (const Translator * translator) : ImplWeightor( translator ) {
}
		       
ImplNoWeightor::~ImplNoWeightor () {
}

ImplNoWeightor::ImplNoWeightor (const ImplNoWeightor & src ) : ImplWeightor( src ) {
}

//--------------------------------------------------------------------------------------------------------------------------------
SequenceWeight * ImplNoWeightor::calculateWeights( const MultipleAlignment & src ) const 
{
  debug_func_cerr(5);


  int nsequences = src.getWidth();

  SequenceWeight * weights = new SequenceWeight[nsequences];
  
  for (int i = 0; i < nsequences; i++) 
    weights[i] = 1;

  return weights;
}

} // namespace alignlib
