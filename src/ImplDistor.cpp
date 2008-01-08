//--------------------------------------------------------------------------------
// Project LibAlign
//
// Copyright (C) 2000 Andreas Heger All rights reserved
//
// Author: Andreas Heger <heger@ebi.ac.uk>
//
// $Id: ImplDistor.cpp,v 1.1.1.1 2002/07/08 21:20:17 heger Exp $
//--------------------------------------------------------------------------------    


#include <iostream>
#include <iomanip>

#include "MultipleAlignment.h"
#include "alignlib.h"
#include "ImplDistor.h"
#include "AlignlibDebug.h"
#include "AlignException.h"
#include "PhyloMatrix.h"

using namespace std;

namespace alignlib {

//---------------------------------------------------------< constructors and destructors >--------------------------------------
ImplDistor::ImplDistor () : mLength(0) {
}
		       
ImplDistor::~ImplDistor () {
}

ImplDistor::ImplDistor (const ImplDistor & src ) : 
    mLength(src.mLength) {
}

//--------------------------------------------------------------------------------------------------------------------------------
PhyloMatrix * ImplDistor::calculateMatrix( PhyloMatrix * matrix, const alignlib::HMultipleAlignment multali) const {

    PhyloMatrixSize i, j;
    
    PhyloMatrixSize width = multali->getWidth();

    if (matrix->getWidth() != width)
	throw AlignException( "Multiple alignment and matrix have different size in ImplDistor::operator()");

    for (i = 0; i < width - 1; i++) 
      for (j = i + 1; j < width; j++) 
	(*matrix)(i, j) = calculateDistance( (*multali)[i], (*multali)[j] );

    return matrix;
} 
    

//--------------------------------------------------------------------------------------------------------------------------------

} // namespace alignlib
