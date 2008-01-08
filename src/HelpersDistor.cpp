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

#include "alignlib.h"
#include "alignlib_fwd.h"
#include "Distor.h"
#include "HelpersDistor.h"

using namespace std;

namespace alignlib 
{
  
  static HDistor DEFAULT_DISTOR(makeDistorClustal());  

  /** gets the default Distor object */
  const HDistor getDefaultDistor() 
  {
    return DEFAULT_DISTOR;
  }
 
  /** sets the default Distor object */
  void setDefaultDistor( const HDistor & distor ) 
  {
    DEFAULT_DISTOR = distor;
  }            

} // namespace alignlib
