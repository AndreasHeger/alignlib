//--------------------------------------------------------------------------------
// Project alignlib
//
// Copyright (C) 2000 Andreas Heger All rights reserved
//
// Author: Andreas Heger <heger@ebi.ac.uk>
//
// $Id: Matrix.cpp,v 1.1.1.1 2002/07/08 21:20:17 heger Exp $
//--------------------------------------------------------------------------------

#include <iostream>
#include "PhyloMatrix.h"

using namespace std;

namespace alignlib {

//---------------------------------------------------------< constructors and destructors >--------------------------------------
PhyloMatrix::PhyloMatrix ()
{
}
		       
PhyloMatrix::~PhyloMatrix ()
{
}

PhyloMatrix::PhyloMatrix (const PhyloMatrix & src )
{
}

std::ostream & operator<<( std::ostream & output, const PhyloMatrix & src)
{
  src.write( output );
  return output;
}
  
std::istream & operator>>( std::istream & input, PhyloMatrix & target)
{
  target.read( input );
  return input;
} 

//---------------------------------------------------------< Input/Output routines >---------------------------------------------
} /* namespace alignlib */
