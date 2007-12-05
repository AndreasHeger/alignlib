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

namespace alignlib {

    class Treetor;
    class Distor;

    /** Helper functions for class Alignata:
	
	1. factory functions
	
	2. accessor functions for default objects
	
	3. convenience functions
    */
    
 /* -------------------------------------------------------------------------------------------------------------------- */
 /* 1. factory functions */

    /** make a Treetor using Distance matrices */
    Treetor * makeTreetorDistanceLinkage( LinkageType method = UPGMA, const Distor * distor = NULL);

    /** make a Treetor using Distance matrices and the neighbour-joining algorithm */
    Treetor * makeTreetorDistanceNJ( const Distor * distor = NULL);    

 /* -------------------------------------------------------------------------------------------------------------------- */
 /* 2. accessor functions for default objects */
    
    /** gets the default Treetor object */
    const Treetor * getDefaultTreetor();
      
    /** sets the default Treetor object
     * 
     * The library obtains ownership of the supplied treetor */
    void setDefaultTreetor( Treetor * treetor );
    
    /* -------------------------------------------------------------------------------------------------------------------- */
    /* 3. convenience functions */


}

#endif	/* HELPERS_TREETOR_H */
