/*
  alignlib - a library for aligning protein sequences

  $Id: ImplDottorSimilarityTuples.cpp,v 1.2 2004/01/07 14:35:35 aheger Exp $

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
#include <iomanip>
#include <map>
#include <vector>
#include <set>
#include <math.h>
#include "alignlib.h"
#include "AlignlibDebug.h"
#include "Alignandum.h"
#include "SubstitutionMatrix.h"
#include "ImplDottorSimilarityTuples.h"
#include "ImplAlignataMatrix.h"

#ifdef WITH_DMALLOC
#include <dmalloc.h>
#endif

using namespace std;

namespace alignlib {

  // map for mapping strings (tuples) to sequenceposotions.
typedef map<string, vector<Position> > TUPLES;

Dottor * makeDottorSimilarityTuples( const SubstitutionMatrix * matrix, int ktuple = 3) {
  return new ImplDottorSimilarityTuples( matrix, ktuple );
}

//---------------------------------------------------------< constructors and destructors >--------------------------------------
ImplDottorSimilarityTuples::ImplDottorSimilarityTuples ( const SubstitutionMatrix * matrix, int ktuple) : ImplDottorSimilarity(matrix), mKtuple(ktuple) {
}
		       
ImplDottorSimilarityTuples::~ImplDottorSimilarityTuples () {
}

ImplDottorSimilarityTuples::ImplDottorSimilarityTuples (const ImplDottorSimilarityTuples & src ) : ImplDottorSimilarity(src), mKtuple( src.mKtuple) {
}

void ImplDottorSimilarityTuples::calculateNewPairs( const Alignandum * row, const Alignandum * col) const 
{
  debug_func_cerr(5);

 
  TUPLES tuples;
 
  // 1. create map of tuples of length ktuple for row_sequence (this is a map of vectors)
  int row_len = row->getLength();

  // get sequence/ consensus-string, etc..
  std::string row_sequence = row->asString();      
 
  for (int xrow = 1; xrow <= (row_len - mKtuple + 1); xrow++)
    tuples[row_sequence.substr(xrow, mKtuple)].push_back(xrow);
 
#ifdef DEBUG
  cout << "Tuples " << endl;
  for(TUPLES::iterator it = tuples.begin(); it != tuples.end(); it++) {
    cout << (*it).first << endl;
    vector<int>& v = tuples[(*it).first];
    copy(v.begin(), v.end(),
         ostream_iterator<int>(cout, " "));
    cout << endl;
  }
#endif

  // 2. look up tuples in col_sequence and add dots to a set, so that they are unique
  int col_len = col->getLength();
  std::string col_sequence = col->asString();
 
  int len = col_len + 1;
 
  set<Position> newdots;                                                              

  // go through object in col
  for (Position xcol = 1; xcol <= (col_len - mKtuple + 1); xcol++) {

    // build tuple
    string s = col_sequence.substr(xcol, mKtuple);

    // look up tuple
    if ( tuples.find(s) != tuples.end() ) {

      // if tuple was there
      vector<Position>& rows = tuples[s];

      for (vector<Position>::iterator it = rows.begin(); it != rows.end(); it++) {

        for (int i = 0; i < mKtuple; i++) {
          int code = (*it + i) * len + xcol + i;
          newdots.insert( code );
        }
      }
    }
  }
 
  FreeMemory();
  AllocateMemory(MAX_NDOTS);
 
  const MATRIXCOLUMN * matrix = (MATRIXCOLUMN*)mMatrix->GetMatrix();
 
  // add dots from set to dots, nice side-effect: the dots are already sorted by row and then by column !!
  for (set<int>::iterator it = newdots.begin(); it != newdots.end(); it++) {
    int xrow = (*it) / len;
    int xcol = (*it) % len;
    AddDot( xrow, xcol, mSubstitutionMatrix->getScore( row->asResidue(xrow), col->asResidue(xcol) ));
  }
 
};


} // namespace alignlib
