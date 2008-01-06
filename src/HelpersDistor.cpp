//--------------------------------------------------------------------------------
// Project LibPhylo
//
// Copyright (C) 2000 Andreas Heger All rights reserved
//
// Author: Andreas Heger <heger@ebi.ac.uk>
//
// $Id: HelpersDistor.cpp,v 1.2 2004/06/02 12:14:34 aheger Exp $
//--------------------------------------------------------------------------------    

#include <iostream>
#include <string>

#ifdef WITH_DMALLOC
#include <dmalloc.h>
#endif

#include "Distor.h"
#include "HelpersDistor.h"

using namespace std;

namespace alignlib 
{
  
  static std::auto_ptr<Distor>DEFAULT_DISTOR(makeDistorClustal());  

  /** gets the default Distor object */
  const Distor * getDefaultDistor() 
  {
    return &*DEFAULT_DISTOR;
  }
 
  /** sets the default Distor object */
  void setDefaultDistor( Distor * distor ) 
  {
    DEFAULT_DISTOR.reset(distor);
  }            

} // namespace alignlib
