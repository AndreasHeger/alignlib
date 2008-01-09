//--------------------------------------------------------------------------------
// Project LibAlign
//
// Copyright (C) 2000 Andreas Heger All rights reserved
//
// Author: Andreas Heger <heger@ebi.ac.uk>
//
// $Id: test_Matrix.cpp,v 1.2 2004/06/02 12:14:35 aheger Exp $
//--------------------------------------------------------------------------------

/** Test Matrix object
 */

#if HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>

#include <time.h> 

#include "alignlib.h"
#include "alignlib_fwd.h"
#include "PhyloMatrix.h"
#include "HelpersPhyloMatrix.h"

using namespace std;

using namespace alignlib;

int main ()
{

	HPhyloMatrix m = makePhyloMatrixSymmetric(10, 0);

	cout << *m << endl;

	// reading and writing
	(*m)(3,5) = 5;
	cout << (*m)(3,5) << endl;

	cout << *m << endl;

	for (int i = 0; i < 10; i++) {
		(*m)(3,i) = i;
		(*m)(5,i) = 100 + i;
	}

	(*m)(0,0) = 0;
	(*m)(3,5) = 1000;
	cout << *m << endl;

	m->swap( 3, 5);

	cout << "After swapping 3 and 5\n" << *m << endl;

	m->shrink();

	cout << "After shrinking\n" << *m << endl;

	return (EXIT_SUCCESS);
}
