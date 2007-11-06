//--------------------------------------------------------------------------------
// Project alignlib
//
// Copyright (C) 2000 Andreas Heger All rights reserved
//
// Author: Andreas Heger <heger@ebi.ac.uk>
//
// $Id: Tree.cpp,v 1.1.1.1 2002/07/08 21:20:17 heger Exp $
//--------------------------------------------------------------------------------

#include "Tree.h"

using namespace std;

namespace alignlib {

//---------------------------------------------------------< constructors and destructors >--------------------------------------
Tree::Tree () {
}
		       
Tree::~Tree () {
}

Tree::Tree (const Tree & src ) {
}

std::ostream & operator<<( std::ostream & output, const Tree & src) {
  src.Write( output );
  return output;
}

  /* 
std::istream & operator>>( std::istream & input, Alignata & target) {
  target.Read( input );
  return input;
} 
  */

//---------------------------------------------------------< Input/Output routines >---------------------------------------------
} /* namespace alignlib */
