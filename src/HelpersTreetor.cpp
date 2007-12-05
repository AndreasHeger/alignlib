//--------------------------------------------------------------------------------
// Project LibPhylo
//
// Copyright (C) 2000 Andreas Heger All rights reserved
//
// Author: Andreas Heger <heger@ebi.ac.uk>
//
// $Id: HelpersTreetor.cpp,v 1.2 2004/06/02 12:14:34 aheger Exp $
//--------------------------------------------------------------------------------    

#include <iostream>
#include <string>

#ifdef WITH_DMALLOC
#include <dmalloc.h>
#endif

#include "Treetor.h"
#include "HelpersTreetor.h"

using namespace std;

namespace alignlib {
  
  static const Treetor * DEFAULT_TREETOR = makeTreetorDistanceLinkage();  

  /** gets the default Treetor object */
  const Treetor * getDefaultTreetor() {
    return DEFAULT_TREETOR;
  }
 
  /** sets the default Treetor object */
  void setDefaultTreetor( Treetor * treetor ) 
  {	
	  if (DEFAULT_TREETOR != NULL)
		 delete DEFAULT_TREETOR;
    DEFAULT_TREETOR = treetor;
  }            

} // namespace alignlib
