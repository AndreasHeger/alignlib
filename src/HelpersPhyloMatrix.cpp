//--------------------------------------------------------------------------------
// Project LibPhylo
//
// Copyright (C) 2000 Andreas Heger All rights reserved
//
// Author: Andreas Heger <heger@ebi.ac.uk>
//
// $Id: HelpersMatrix.cpp,v 1.2 2004/06/02 12:14:34 aheger Exp $
//--------------------------------------------------------------------------------    

#include <iostream>
#include <string>

#include "alignlib.h"
#include "alignlib_fwd.h"
#include "AlignlibDebug.h"

#include "HelpersPhyloMatrix.h"
#include "PhyloMatrix.h"

using namespace std;

namespace alignlib 
{

  /** copy the contents of source element wise into the PhyloMatrix */
  HPhyloMatrix & fillPhyloMatrix( HPhyloMatrix & dest, 
		  PhyloMatrixValue * source) 
  {
	  debug_func_cerr( 5 );
    
      PhyloMatrixSize row, col, index = 0;
    
      for (row = 0; row < dest->getWidth(); row++)
    	  for (col = 0; col < dest->getWidth(); col++) 
    	  {
    		  cout << (*dest)(0,1);
    		  (*dest)(row,col) = source[index++];
    		  cout << " " << row << " " << col << (*dest)(0,1) << endl;
    	  }

      cout << *dest << endl;
      return dest;
  }
  
} // namespace alignlib
