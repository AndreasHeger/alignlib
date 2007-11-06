//--------------------------------------------------------------------------------
// Project LibPhylo
//
// Copyright (C) 2000 Andreas Heger All rights reserved
//
// Author: Andreas Heger <heger@ebi.ac.uk>
//
// $Id: HelpersTree.h,v 1.1.1.1 2002/07/08 21:20:17 heger Exp $
//--------------------------------------------------------------------------------


#if HAVE_CONFIG_H
#include <config.h>
#endif

#ifndef HELPERS_TREE_H
#define HELPERS_TREE_H 1

#include <iosfwd>
#include <string>
#include "alignlib.h"

namespace alignlib {

    class Tree;
    class Distor;

    /** Helper functions for class Alignata:
	
	1. factory functions
	
	2. accessor functions for default objects
	
	3. convenience functions
    */
    
    /* -------------------------------------------------------------------------------------------------------------------- */
    /* 1. factory functions */

    /** create an empty tree */
    Tree * makeTree( Node num_leaves = 0);

    /* -------------------------------------------------------------------------------------------------------------------- */
    /* 2. accessor functions for default objects */
    
    /* -------------------------------------------------------------------------------------------------------------------- */
    /* 3. convenience functions */

    /** print a tree in NewHampshire format */
    void writeNewHampshire( std::ostream & output, const Tree * tree, const std::vector<std::string> * labels = NULL );
 
}

#endif	/* HELPERS_TREE_H */
