//--------------------------------------------------------------------------------
// Project LibPhylo
//
// Copyright (C) 2000 Andreas Heger All rights reserved
//
// Author: Andreas Heger <heger@ebi.ac.uk>
//
// $Id: HelpersTreetor.h,v 1.1.1.1 2002/07/08 21:20:17 heger Exp $
//--------------------------------------------------------------------------------


#if HAVE_CONFIG_H
#include <config.h>
#endif

#ifndef HELPERS_TREETOR_H
#define HELPERS_TREETOR_H 1

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

    /** make a Treetor using Distance matrices */
	HTreetor makeTreetorDistanceLinkage(
		const HDistor & distor,     		
		LinkageType method = UPGMA ); 

    /** make a Treetor using Distance matrices and the neighbour-joining algorithm */
    HTreetor makeTreetorDistanceNJ( const HDistor & distor );    

    /* -------------------------------------------------------------------------------------------------------------------- */
    /* 2. accessor functions for default objects */

    DEFINE_DEFAULT( HTreetor, 
    		getDefaultTreetor,
    		setDefaultTreetor );
    
    /* -------------------------------------------------------------------------------------------------------------------- */
    /* 3. convenience functions */
}

#endif	/* HELPERS_TREETOR_H */
