/*
  alignlib - a library for aligning protein sequences

  $Id: ImplRegularizorDirichletFade.cpp,v 1.2 2004/01/07 14:35:35 aheger Exp $

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
#include "ImplRegularizorDirichletFade.h"
#include <math.h>
#include <iostream>

namespace alignlib {

    /** factory functions */
Regularizor * makeRegularizorDirichletFade() { return new ImplRegularizorDirichletFade();}

/* This cutoff determines the number of counts that are needed to switch from
   Dirichlet-mixtures to 
*/
#define FADE_CUTOFF 10
//---------------------------------------------------------< constructors and destructors >--------------------------------------
ImplRegularizorDirichletFade::ImplRegularizorDirichletFade () : ImplRegularizorDirichlet() {
}
  
//--------------------------------------------------------------------------------------------------------------------------------
ImplRegularizorDirichletFade::~ImplRegularizorDirichletFade () {
}

//-------------------------------------------------------------------------------------------------------------------------------
ImplRegularizorDirichletFade::ImplRegularizorDirichletFade (const ImplRegularizorDirichletFade & src ) : ImplRegularizorDirichlet( src ) {
}

//-------------------------------------------------------------------------------------------------------
/** see Kimmen's PhD-Thesis, pp42ff.
      I use equations 6.38 and 6.39 for calculating the pi
      using the logarithm of the beta-function, as Kimmen suggests.
      This method calculates the Xi and stores it in the profile. You have
      to call NormalizeColumn to get the correct estimated probabilites.
*/      
void ImplRegularizorDirichletFade::fillFrequencies( FrequencyColumn * frequencies, 
						    const CountColumn * counts, 
						    const Position length ) const {
    Position column;
    Count ntotal;

    int i,j;
    double Xi[PROFILEWIDTH];
    
    // helper variables for the conversion of log to floating point
    TYPE_BETA_DIFFERENCES beta_differences;
 
    for (column = 0; column < length; column++) {

	const Count* n = counts[column];
	
	ntotal = 0;

	// get ntotal = number of observations
	for (i = 0; i < PROFILEWIDTH; i++)
	    ntotal += n[i];
	
	// if the number of total observation is less than the fade cutoff, do the normal calculation. I have pasted and copied
	// it, maybe externalize it sometime later..

	if (ntotal < FADE_CUTOFF) {
	    // use normal Dirichlet-Mixture
	    fillColumn( beta_differences, n. ntotal );

	}  else {
	    // calculate raw frequencies:
	    for (i = 0; i < PROFILEWIDTH; i++)
		frequencies[column][i] = (Frequency)(n[i] / ntotal);
	    
	}
    }
}

} // namespace alignlib
