//--------------------------------------------------------------------------------
// Project LibAlign
//
// Copyright (C) 2000 Andreas Heger All rights reserved
//
// Author: Andreas Heger <heger@ebi.ac.uk>
//
// $Id: ImplDistorDummy.cpp,v 1.1.1.1 2002/07/08 21:20:17 heger Exp $
//--------------------------------------------------------------------------------    


#include <iostream>
#include "ImplDistorDummy.h"
#include "PhyloMatrix.h"
#include "AlignlibDebug.h"

using namespace std;

namespace alignlib {

//-------------------------> factory functions <-------------------------------------------------------------------------------
Distor * makeDistorDummy( const PhyloMatrix * matrix) {
  return new ImplDistorDummy( matrix);
}

//---------------------------------------------------------< constructors and destructors >--------------------------------------
ImplDistorDummy::ImplDistorDummy ( const PhyloMatrix * matrix) : ImplDistor(), mMatrix(matrix) {
}
		       
ImplDistorDummy::~ImplDistorDummy () 
{
	debug_func_cerr( 5 );
}

ImplDistorDummy::ImplDistorDummy (const ImplDistorDummy & src ) : ImplDistor( src ), mMatrix(src.mMatrix) {
}


//--------------------------------------------------------------------------------------------------------------------------------
PhyloMatrixValue ImplDistorDummy::getMaximumPossibleDistance() const 
{
    return mMatrix->getMaximum();
}

//--------------------------------------------------------------------------------------------------------------------------------
PhyloMatrix * ImplDistorDummy::calculateMatrix( PhyloMatrix * matrix, const alignlib::MultipleAlignment * multali) const 
{
	debug_func_cerr( 5 );

    PhyloMatrixSize i, j;
    PhyloMatrixSize width = mMatrix->getWidth();
    
    matrix->setWidth( width );

    // This is definitely not efficient
    for (i = 0; i < width - 1; i++) 
	for (j = i + 1; j < width; j++) 
	    (*matrix)(i, j) = (*mMatrix)(i, j);
    
    return matrix;
} 
    
//--------------------------------------------------------------------------------------------------------------------------------
PhyloMatrixValue ImplDistorDummy::calculateDistance( const std::string & s_row_1, const std::string & s_row_2) const 
{
	debug_func_cerr( 5 );
  
  return 0;
}


} // namespace alignlib
