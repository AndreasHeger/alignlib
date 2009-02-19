//--------------------------------------------------------------------------------
// Project LibAlign
//
// Copyright (C) 2000 Andreas Heger All rights reserved
//
// Author: Andreas Heger <heger@ebi.ac.uk>
//
// $Id: test_Tree.cpp,v 1.2 2004/06/02 12:14:35 aheger Exp $
//--------------------------------------------------------------------------------

/** Test alignata objects
 */

#if HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>
#include <iterator>
#include <fstream>
#include <algorithm>
#include <vector>
#include <cassert>
#include <time.h>

#include "alignlib.h"
#include "alignlib_fwd.h"
#include "HelpersTree.h"
#include "Tree.h"

using namespace std;
using namespace alignlib;

int main ()
{

	HTree tree(makeTree(4));

	cout << *tree << endl;

	cout << "Joining nodes 0 and 1 gives " << tree->joinNodes( 0, 1, 1, 2) << endl;
	cout << "Joining nodes 2 and 3 gives " << tree->joinNodes( 2, 3, 1, 2) << endl;
	cout << "Joining nodes 4 and 5 gives " << tree->joinNodes( 4, 5, 1, 2) << endl;

	cout << *tree << endl;

	writeNewHampshire( cout, tree, HStringVector( new StringVector()) );

	cout << "Output from setRoot:" << tree->setRoot( 0, 4, 1) << endl;

	cout << *tree << endl;

	// writeNewHampshire( cout, tree );

	HNodeVector v(tree->getNodesLeaves());
	cout << "Leaves-traversal:" << endl;
	std::copy( v->begin(), v->end(), std::ostream_iterator< Node >( std::cout, " " ));
	cout << endl;

	v = tree->getNodesDepthFirstVisit();
	cout << "DFS-traversal visit:" << endl;
	std::copy( v->begin(), v->end(), std::ostream_iterator< Node >( std::cout, " " ));
	cout << endl;

	v = tree->getNodesDepthFirstFinish();
	cout << "DFS-traversal finish:" << endl;
	std::copy( v->begin(), v->end(), std::ostream_iterator< Node >( std::cout, " " ));
	cout << endl;

	v = tree->getNodesBreadthFirstVisit();
	cout << "BFS-traversal visit:" << endl;
	std::copy( v->begin(), v->end(), std::ostream_iterator< Node >( std::cout, " " ));
	cout << endl;

	{
		std::cout << "copying" << std::endl;
		HTree t2 = tree->getClone();
		assert( t2->getRoot() == tree->getRoot() );
		assert( t2->getNumLeaves() == tree->getNumLeaves() );
	}

	return (EXIT_SUCCESS);
}
