/*
  alignlib - a library for aligning protein sequences

  $Id: ImplSubstitutionMatrixAA.cpp,v 1.3 2004/06/02 11:52:49 aheger Exp $

  Copyright (C) 2004 Andreas Heger
  
  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License
  as published by the Free Software Foundation; either version 2
  of the License, or (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
*/


#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdlib.h>
#include "alignlib.h"
#include "AlignlibDebug.h"
#include "AlignException.h"
#include "ImplSubstitutionMatrixAA.h"

#ifdef WITH_DMALLOC
#include <dmalloc.h>
#endif

using namespace std;

namespace alignlib {


//------------------------------------------------------------------------------------------------------------
/** defintion of global object */

//------------> factory functions <------------------------------------------------------------------------------------------------
SubstitutionMatrix * makeSubstitutionMatrixAA( ScoreColumn * matrix, bool this_own) {
  return new ImplSubstitutionMatrixAA( matrix, this_own );
}

//----------------------------------------------------------------------------------------------------------------------------
SubstitutionMatrix * readSubstitutionMatrixAA( const char * filename ) {
    // simple reading routine for sustitution matrix
  
  ifstream fin( filename);  
  if (!fin)       
    throw AlignException("Could not open file in readSubstitutionMatrixAA");
  
  ScoreColumn * matrix = new ScoreColumn[MATRIXWIDTH_AA];
  
  Score * p = (Score*)matrix;
  int max_iterations = MATRIXWIDTH_AA * MATRIXWIDTH_AA;
  
  while (max_iterations > 0 && !fin.eof()) {                          // while (fin) does not work!!       
    fin >> *p;
    p++;
    max_iterations--;
  }
    
  if (max_iterations != 0 )
    throw AlignException("Read incomplete matrixin readSubstitionMatrixAA");

  fin.close();  
  
  return new ImplSubstitutionMatrixAA( matrix, true );
}


//---------------------------------------------------------< constructors and destructors >--------------------------------------
ImplSubstitutionMatrixAA::ImplSubstitutionMatrixAA( ScoreColumn * matrix, bool this_own ) :
    ImplSubstitutionMatrix( (Score*)matrix, MATRIXWIDTH_AA, MATRIXWIDTH_AA, this_own) {

    mMatrixAA = matrix;
    mData.mNumRows = MATRIXWIDTH_AA;
    mData.mNumCols = MATRIXWIDTH_AA;
    mData.mAAMatrixPointer = matrix;
}

ImplSubstitutionMatrixAA::~ImplSubstitutionMatrixAA () {
}
  
ImplSubstitutionMatrixAA::ImplSubstitutionMatrixAA (const ImplSubstitutionMatrixAA & src ) : ImplSubstitutionMatrix( src) {
}

//--------------------------------------------------------------------------------------------------------------------------------
const ImplSubstitutionMatrixAAData & ImplSubstitutionMatrixAA::getData() const { return mData;}

//--------------------------------------------------------------------------------------------------------------------------------
Score ImplSubstitutionMatrixAA::getScore( Residue row, Residue col) const {
    return mMatrixAA[row][col];
}


} // namespace alignlib
