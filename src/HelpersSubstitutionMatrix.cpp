/*
  alignlib - a library for aligning protein sequences

  $Id: HelpersSubstitutionMatrix.cpp,v 1.2 2004/01/07 14:35:33 aheger Exp $

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

#include "SubstitutionMatrix.h"
#include "HelpersSubstitutionMatrix.h"

using namespace std;

namespace alignlib 
{

    static ScoreColumn blosum62[MATRIXWIDTH_AA] = {
/*   A   C   D   E   F   G   H   I   K   L   M   N   P   Q   R   S   T   V   W   Y   X*/
{    4,  0, -2, -1, -2,  0, -2, -1, -1, -1, -1, -2, -1, -1, -1,  1,  0,  0, -3, -2, -4,}, /* A*/
{    0,  9, -3, -4, -2, -3, -3, -1, -3, -1, -1, -3, -3, -3, -3, -1, -1, -1, -2, -2, -4,}, /* C*/
{   -2, -3,  6,  2, -3, -1, -1, -3, -1, -4, -3,  1, -1,  0, -2,  0, -1, -3, -4, -3, -4,}, /* D*/
{   -1, -4,  2,  5, -3, -2,  0, -3,  1, -3, -2,  0, -1,  2,  0,  0, -1, -2, -3, -2, -4,}, /* E*/
{   -2, -2, -3, -3,  6, -3, -1,  0, -3,  0,  0, -3, -4, -3, -3, -2, -2, -1,  1,  3, -4,}, /* F*/
{    0, -3, -1, -2, -3,  6, -2, -4, -2, -4, -3,  0, -2, -2, -2,  0, -2, -3, -2, -3, -4,}, /* G*/
{   -2, -3, -1,  0, -1, -2,  8, -3, -1, -3, -2,  1, -2,  0,  0, -1, -2, -3, -2,  2, -4,}, /* H*/
{   -1, -1, -3, -3,  0, -4, -3,  4, -3,  2,  1, -3, -3, -3, -3, -2, -1,  3, -3, -1, -4,}, /* I*/
{   -1, -3, -1,  1, -3, -2, -1, -3,  5, -2, -1,  0, -1,  1,  2,  0, -1, -2, -3, -2, -4,}, /* K*/
{   -1, -1, -4, -3,  0, -4, -3,  2, -2,  4,  2, -3, -3, -2, -2, -2, -1,  1, -2, -1, -4,}, /* L*/
{   -1, -1, -3, -2,  0, -3, -2,  1, -1,  2,  5, -2, -2,  0, -1, -1, -1,  1, -1, -1, -4,}, /* M*/
{   -2, -3,  1,  0, -3,  0,  1, -3,  0, -3, -2,  6, -2,  0,  0,  1,  0, -3, -4, -2, -4,}, /* N*/
{   -1, -3, -1, -1, -4, -2, -2, -3, -1, -3, -2, -2,  7, -1, -2, -1, -1, -2, -4, -3, -4,}, /* P*/
{   -1, -3,  0,  2, -3, -2,  0, -3,  1, -2,  0,  0, -1,  5,  1,  0, -1, -2, -2, -1, -4,}, /* Q*/
{   -1, -3, -2,  0, -3, -2,  0, -3,  2, -2, -1,  0, -2,  1,  5, -1, -1, -3, -3, -2, -4,}, /* R*/
{    1, -1,  0,  0, -2,  0, -1, -2,  0, -2, -1,  1, -1,  0, -1,  4,  1, -2, -3, -2, -4,}, /* S*/
{    0, -1, -1, -1, -2, -2, -2, -1, -1, -1, -1,  0, -1, -1, -1,  1,  5,  0, -2, -2, -4,}, /* T*/
{    0, -1, -3, -2, -1, -3, -3,  3, -2,  1,  1, -3, -2, -2, -3, -2,  0,  4, -3, -1, -4,}, /* V*/
{   -3, -2, -4, -3,  1, -2, -2, -3, -3, -2, -1, -4, -4, -2, -3, -3, -2, -3, 11,  2, -4,}, /* W*/
{   -2, -2, -3, -2,  3, -3,  2, -1, -2, -1, -1, -2, -3, -1, -2, -2, -2, -1,  2,  7, -4,}, /* Y*/
{   -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4,}, /* X */
    };

  static SubstitutionMatrix * DEFAULT_SUBSTITUTIONMATRIX = makeSubstitutionMatrixAA( blosum62 );
  
  /** gets the default SubstitutionMatrix object */ 
  const SubstitutionMatrix * getDefaultSubstitutionMatrix() 
  {
    return DEFAULT_SUBSTITUTIONMATRIX;
  }

  /** sets the default SubstitutionMatrix object */
  void setDefaultSubstitutionMatrix( SubstitutionMatrix * matrix ) 
  {
	  if (DEFAULT_SUBSTITUTIONMATRIX != NULL)
		  delete DEFAULT_SUBSTITUTIONMATRIX;
	  DEFAULT_SUBSTITUTIONMATRIX = matrix;
  }
 
/** create the identity substitution matrix. Identities score as 1, mismatches as -1.
 */

SubstitutionMatrix * makeSubstitutionMatrixAAIdentity( Score match,
						       Score mismatch ) {
  ScoreColumn * matrix = new ScoreColumn[MATRIXWIDTH_AA];
  
  unsigned int row, col;
  for (row = 0; row < MATRIXWIDTH_AA; row++) {
    for (col = 0; col < MATRIXWIDTH_AA; col++) 
      matrix[row][col] = mismatch;
    matrix[row][row] = match;
  }
  
  return makeSubstitutionMatrixAA( matrix, true );

}

} // namespace alignlib









