/*
  alignlib - a library for aligning protein sequences

  $Id: ImplNoRegularizor.cpp,v 1.2 2004/01/07 14:35:35 aheger Exp $

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
#include "Regularizor.h"
#include "HelpersRegularizor.h"
#include "ImplNoRegularizor.h"

#ifdef WITH_DMALLOC
#include <dmalloc.h>
#endif

using namespace std;

namespace alignlib {

  /** factory functions */
  Regularizor * makeNoRegularizor() { return new ImplNoRegularizor();}


//---------------------------------------------------------< constructors and destructors >--------------------------------------
ImplNoRegularizor::ImplNoRegularizor () {
}
		       
ImplNoRegularizor::~ImplNoRegularizor () {
}

ImplNoRegularizor::ImplNoRegularizor (const ImplNoRegularizor & src ) {
}

//--------------------------------------------------------------------------------------------------------------------------------
void ImplNoRegularizor::fillFrequencies( FrequencyColumn * frequencies, 
					 const CountColumn * counts, 
					 const Position length ) const 
{
  debug_func_cerr(5);


  // simply calculate frequencies
 
  Position column;
  Count ntotal;
  int i;

  for (column = 1; column <= length; column++) {
    ntotal = 0;
 
    for (i = 0; i < PROFILEWIDTH; i ++)
      ntotal += counts[column][i];

    if (ntotal == 0)
      ntotal = 1;
    
    for (i = 0; i < PROFILEWIDTH; i++) {
      Frequency f = (Frequency)((Frequency)counts[column][i] / (Frequency)ntotal);
      frequencies[column][i] = f;
    }
    
  }                    

}

} // namespace alignlib
