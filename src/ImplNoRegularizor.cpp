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

using namespace std;

namespace alignlib 
{

/** factory functions */
Regularizor * makeNoRegularizor() 
{ 
	return new ImplNoRegularizor();
}


//---------------------------------------------------------< constructors and destructors >--------------------------------------
ImplNoRegularizor::ImplNoRegularizor () {
}
		       
ImplNoRegularizor::~ImplNoRegularizor () {
}

ImplNoRegularizor::ImplNoRegularizor (const ImplNoRegularizor & src ) 
{
}

//--------------------------------------------------------------------------------------------------------------------------------
void ImplNoRegularizor::fillFrequencies( FrequencyMatrix * frequencies, 
					 const CountMatrix * counts ) const
					 {

  debug_func_cerr(5);
  
  // simply calculate frequencies
 
  Position column;
  Count ntotal;
  int i;

  Position width = frequencies->getNumCols();
  Position length = frequencies->getNumRows();
  
  for (column = 0; column < length; ++column) 
  {
    ntotal = 0;
 
    const Count * counts_column = counts->getRow(column);
    Frequency * frequency_column = frequencies->getRow(column);
    
    for (i = 0; i < width; i ++)
      ntotal += counts_column[i];

    if (ntotal == 0)
      ntotal = 1;
    
    for (i = 0; i < width; i++) 
    {
      Frequency f = (Frequency)((Frequency)counts_column[i] / (Frequency)ntotal);
      frequency_column[i] = f;
    }
    
  }                    

}

} // namespace alignlib
