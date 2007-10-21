/*
  alignlib - a library for aligning protein sequences

  $Id: ImplWeightor.cpp,v 1.2 2004/01/07 14:35:36 aheger Exp $

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
#include "ImplWeightor.h"
#include "HelpersWeightor.h"

#ifdef WITH_DMALLOC
#include <dmalloc.h>
#endif

using namespace std;

namespace alignlib {

#define MIN_WEIGHT 0.0001

//---------------------------------------------------------< constructors and destructors >--------------------------------------
ImplWeightor::ImplWeightor ( const Translator * translator) : mTranslator(translator){
}
		       
ImplWeightor::ImplWeightor (const ImplWeightor & src ) : Weightor(src), mTranslator( src.mTranslator ) {
}

ImplWeightor::~ImplWeightor () 
{
  debug_func_cerr(5);

}

//--------------------------------------------------------------------------------------------------------------------------------
void ImplWeightor::rescaleWeights( SequenceWeights * weights, int nsequences, SequenceWeight value) const 
{
  debug_func_cerr(5);


    if (value == 0)
      value = nsequences;

    //--------------> rescale weights
    double total = 0;
    int i;
    SequenceWeights & w = *weights;
    
    for ( i = 0; i < nsequences; i++) {
	if (w[i] < MIN_WEIGHT) 
	    w[i] = MIN_WEIGHT;               //!! to do: some warnings, exception handling, ...
	total += w[i];
    }
    
    SequenceWeight factor = value / total;
 
    for ( i = 0; i < nsequences; i++) w[i] *= factor;
}


} // namespace alignlib
