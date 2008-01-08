//--------------------------------------------------------------------------------
// Project LibPhylo
//
// Copyright (C) 2000 Andreas Heger All rights reserved
//
// Author: Andreas Heger <heger@ebi.ac.uk>
//
// $Id: HelpersDistor.h,v 1.1.1.1 2002/07/08 21:20:17 heger Exp $
//--------------------------------------------------------------------------------


#if HAVE_CONFIG_H
#include <config.h>
#endif

#ifndef HELPERS_DISTOR_H
#define HELPERS_DISTOR_H 1

#include <iosfwd>
#include <string>
#include "alignlib.h"
#include "alignlib_fwd.h"
#include "alignlib_default.h"

namespace alignlib 
{
    /** Helper functions for class Alignment:
	
	1. factory functions
	
	2. accessor functions for default objects
	
	3. convenience functions
    */
    
 /* -------------------------------------------------------------------------------------------------------------------- */
 /* 1. factory functions */
  
    /** return an Distor object, that calculates distances between sequences in a multiple alignment based on Kimura */
    HDistor makeDistorKimura();

    /** return an Distor object, that calculates distances between sequences in a multiple alignment based on Kimura */
    HDistor makeDistorClustal();
    
    /** return a Distor object, that instead of calculating a matrix, copies its own copy of a matrix into a matrix */
    HDistor makeDistorDummy( const HPhyloMatrix & matrix);
  
 /* -------------------------------------------------------------------------------------------------------------------- */
 /* 2. accessor functions for default objects */
    
    DEFINE_DEFAULT( HDistor, getDefaultDistor, setDefaultDistor );
    /* -------------------------------------------------------------------------------------------------------------------- */
    /* 3. convenience functions */


}

#endif	/* HELPERS_DISTOR_H */
