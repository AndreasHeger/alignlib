/*
  alignlib - a library for aligning protein sequences

  $Id: ImplAlignatorFullDP.cpp,v 1.2 2004/01/07 14:35:34 aheger Exp $

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
#include "alignlib.h"
#include "AlignlibDebug.h"
#include "AlignException.h"

#include "SubstitutionMatrix.h"
#include "HelpersSubstitutionMatrix.h"
#include "ImplSubstitutionMatrixAA.h"

#include "Alignata.h"
#include "HelpersAlignata.h"

#include "Alignandum.h"
#include "ImplAlignatorFullDP.h"
#include "Alignator.h"

//implementation of Alignandum-objects
#include "ImplSequence.h"
#include "ImplProfile.h"

#include "ImplTranslator.h" // for direct access to mask_code


#ifdef WITH_DMALLOC
#include <dmalloc.h>
#endif

using namespace std;

namespace alignlib {
  
  /* How to write a fast algorithm:
     My design objective here was to not duplicate the algorithmic code without penalizing too much
     for function indirection. It the way I do it below, there is one indirection for every call to
     calculate a match function, except for sequence-sequence comparisons, where there would have been
     two.
    */
       
//<--------------< Global functions and pointers for fast determination of match score. >-------------------------------------------
// I don't know how to use member functions as function pointers. After all, this is what inheritence is for. A possibility would be
// to automatically subclass an alignator-object.However, I do not like this idea, since this assumes that the parent has information 
// about the child. On the other hand, via the inlining mechanism a function call could be saved. The problem is when you want to change
// the algorithm by overloading align. Then the parent functions () do not know, what the child is. Therefore I use the static 
// functions. Maybe it is possible to separate the algorithm and the type-decision into different classes that interact.
//
// The danger of static functions is that the global pointers are unsafe, i.e. there exist just one copy for all alignator-objects, and
// this code will never be threadsafe.
//-----------------------------------------------------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------------------------------------------------------
ImplAlignatorFullDP::ImplAlignatorFullDP( const SubstitutionMatrix * subst_matrix,
					  Score row_gop, Score row_gep, 
					  Score col_gop, Score col_gep ) :
    ImplAlignator( subst_matrix, row_gop, row_gep, col_gop, col_gep), 
    mTrace(NULL), 
    mRowLast(0), mColLast(0),
    mCC(NULL), 
    mDD(NULL) {
}
  
//----------------------------------------------------------------------------------------------------------------------------------------
ImplAlignatorFullDP::ImplAlignatorFullDP( const ImplAlignatorFullDP & src ) : ImplAlignator( src ) 
{
  debug_func_cerr(5);


    mTrace = NULL;
    mRowLast = 0;
    mColLast = 0;

    mCC    = NULL;
    mDD    = NULL;
    
}

//----------------------------------------------------------------------------------------------------------------------------------------
ImplAlignatorFullDP::~ImplAlignatorFullDP() 
{
  debug_func_cerr(5);


}

//----------------------------------------------------------------------------------------------------------------------------------------
Alignata * ImplAlignatorFullDP::align(const Alignandum * row, const Alignandum * col, Alignata * result) 
{
  debug_func_cerr(5);


  /* try casting down the hiearchy of Alignandum-objects and then using a switch statement 
     register the correct match function. There are different implementations for this:
     1. multiple dispatch, but then everytime I write code for aligning/creating dots, ...
     I have to add new corresponding functions in Alignandum objects. I rather like the code
     in one file, even if it is large.
     2. using an emulated virtual table. This is an elegant solution (see Meyers), but needs
     a lot of coding.
     I do not expect, that there will be many more other Alignandum-objects, so I use the
     dynamics_cast-way.
  */
  startUp(row, col, result);

  performAlignment(row, col, result);
    
  traceBack(row, col, result);
  
  cleanUp(row, col, result);
    
  return result;
}

//----------------------------------------------------------------------------------------------------------------------------------------
Score ImplAlignatorFullDP::getAlignmentScore(const Alignandum * row, const Alignandum * col) 
{
  debug_func_cerr(5);


  Alignata * dummy_alignment = makeAlignataSet();

  startUp(row, col, dummy_alignment);

  performAlignmentWithoutTraceBack(row, col);
    
  cleanUp(row, col, dummy_alignment);
    
  delete dummy_alignment;

  return mScore;
}

//----------------------------------------------------------------------------------------------------------------------------------------
void ImplAlignatorFullDP::startUp(const Alignandum * row, const Alignandum *col, Alignata * ali) {
    ImplAlignator::startUp(row, col, ali);  

    mRowLengthPlus1 = getRowLength() + 1;
    //----------------------------------------------------------------------------------------
    // setup pointers to data locations
    // const SequenceTCO & t1 = dynamic_cast<const SequenceTCO&>(row);
    
    /* try casting down the hiearchy of Alignandum-objects and then using a switch statement 
       register the correct match function. There are different implementations for this:
       1. multiple dispatch, but then everytime I write code for aligning/creating dots, ...
       I have to add new corresponding functions in Alignandum objects. I rather like the code
       in one file, even if it is large.
       2. using an emulated virtual table. This is an elegant solution (see Meyers), but needs
       a lot of coding.
       I do not expect, that there will be many more other Alignandum-objects, so I use the
       dynamics_cast-way.
    */
    const ImplSequence * s1 = dynamic_cast<const ImplSequence*>(row);
    const ImplProfile * p1 = dynamic_cast<const ImplProfile*>(row);

    // const SequenceTCO & t2 = dynamic_cast<const SequenceTCO&>(col);
    const ImplSequence * s2 = dynamic_cast<const ImplSequence*>(col);
    const ImplProfile * p2 = dynamic_cast<const ImplProfile*>(col);

    // setup static pointers to the data locations
    if (s1) {
	mRowSequence = s1->getData().mSequencePointer;
    } else if (p1) {
	mRowProfile     = p1->getData().mProfilePointer;
	mRowFrequencies = p1->getData().mFrequenciesPointer;
    }
    //!! to do: add TCO

    if (s2) {
	mColSequence = s2->getData().mSequencePointer;
    } else if (p2) {
	mColProfile     = p2->getData().mProfilePointer;
	mColFrequencies = p2->getData().mFrequenciesPointer;
    }

    //!! to do: add TCO
    if (s1 && s2) {
	mMatchFunction = &ImplAlignatorFullDP::matchSequenceSequence;
    }

    if (p1 && p2) {
	mMatchFunction = &ImplAlignatorFullDP::matchProfileProfile;
    }

    if (s1 && p2) {
	mMatchFunction = &ImplAlignatorFullDP::matchSequenceProfile;
    }
  
    if (p1 && s2) {
	mMatchFunction = &ImplAlignatorFullDP::matchProfileSequence;
    }

    //----------------------------------------------------------------------------------------
    //--------------------------------> setup mTraceback-matrix
    int tb_size = (getRowLength() + 1) * ( getColLength() + 1);

    mTrace = new int[tb_size];

    for (int i = 0; i < tb_size; i++)
	mTrace[i] = TB_STOP;
    
    mCC    = new Score[getColLength() + 1];
    mDD    = new Score[getColLength() + 1];
    
    mRowLast = 0;
    mColLast = 0;
    mScore = 0;

#ifdef DEBUG
    // cout << "Substitution-Matrix used" << endl << *mSubstitutionMatrix;
#endif
}
//----------------------------------------------------------------------------------------------------------------------------------------
void ImplAlignatorFullDP::cleanUp(const Alignandum * row, const Alignandum *col, Alignata * ali) {
    
    if (mCC != NULL)    { delete [] mCC; mCC = NULL; }
    if (mDD != NULL)    { delete [] mDD; mDD = NULL; }
    if (mTrace != NULL) { delete [] mTrace; mTrace = NULL; }
    
    ImplAlignator::cleanUp(row, col, ali);

}       

//-------------------------------------< BackTrace >----------------------------------------------------------------------
// wrapping around for col but not for row, because otherwise there could be an infinite loop.
#define PREVCOL { if (--col < 1) col = mColLength; }  
void ImplAlignatorFullDP::traceBack( const Alignandum * prow, const Alignandum * pcol, Alignata * result) 
{
  debug_func_cerr(5);

 
    int row, col;
    int t;
    
    Position col_length = getColLength();
#ifdef DEBUG
    Position row_length = getRowLength();
    cout << "Trace matrix" << endl;
    cout << setw(8) << "";
    for (col = 0; col <= col_length; col++) cout << setw(4) << col;
    cout << endl;
    for (row = 0; row <= row_length; row++) {
      cout << setw(8) << row;
      for (col = 0; col <= col_length; col++) {
	cout << setw(4);
	switch (mTrace[getTraceIndex(row,col)]) {
	case TB_STOP:      cout << "o"; break;
	case TB_DELETION:  cout << "<" ; break;
	case TB_INSERTION: cout << "^" ; break;
	case TB_MATCH:     cout << "=" ; break;
	case TB_WRAP:	   cout << "@" ; break;
	default: cout << "#"; break;
	}
      }
      cout << endl;
    }
    cout << "-----------------------" << endl;
    cout << " mRowLast " << mRowLast << "mColLast " << mColLast << endl;
#endif
 
    col = mColLast;
    row = mRowLast;
    int ngaps = 0;

    t = mTrace[getTraceIndex(row,col)];

    while ( t != TB_STOP && row > 0) 
{
  debug_func_cerr(5);

        switch (t) {
	case TB_DELETION  :
	    col--;
	    ngaps++;
	    if (col < 1) 
		row--;
	    break;
        case TB_INSERTION :
	    ngaps++;
	    row--;
	    break;
        case TB_MATCH     :
	  result->addPair( new ResiduePAIR( row, col,(this->*mMatchFunction)( row, col)));
	    row--;
	    col--;
	  break;
	case TB_WRAP :
	  col = col_length;
	  break;
        default:
          throw AlignException("Unknown matrix command in TraceBack");
          break;
        }
        t = mTrace[getTraceIndex(row,col)];
    }
    result->setScore ( mScore );
}  

//---------------------------------< the actual alignment algorithm >-------------------------------------------
void ImplAlignatorFullDP::performAlignment( const Alignandum * prow, const Alignandum * pcol, Alignata * ali) {

  register int   row, col;

  Position row_length = getRowLength();
  Position col_length = getColLength();

  Score row_gop = getRowGop();
  Score row_gep = getRowGep();
  Score col_gop = getColGop();
  Score col_gep = getColGep();

  mRowLast = 0;
  mColLast = 0;
  mScore = 0;

  Score c, e, d, s;                  // helper variables
  
  Score row_m = row_gop + row_gep;
  Score col_m = col_gop + col_gep;
  
  //----------------------------> Initialise affine penalty arrays <-------------------------------
  mCC[0] = 0;
  for (col = 1; col <= col_length; col++)
    {
      mCC[col]   = 0;
      mDD[col]   = row_gop;                               // score for horizontal gap opening
    }
 
  mCC[col_length]   = col_gop;
    
  //----------------------------> Calculate dynamic programming matrix <----------------------------
  //----------------------------> iterate over rowumns <--------------------------------------------
  for (row = 1; row <= row_length; row++)
    {
      
      // this part is different from the ordinary alignment
      s = mCC[0];
      mCC[0] = c = 0;
      
      e = col_gop;                                        // penalty for opening a vertical gap
	
      //-------------------------> iterate over cols <------------------------------------------------
      for (col = 1; col <= col_length; col++)
	{
	  
	  //---------------------------> calculate scores <--------------------------------------------
	  // c contains score of cell above
	  // s contains score for cell [row, col-1]
	  // e is better of: score for opening a vertical gap or score for extending a vertical gap: use col-gap-penalties
	  if ((c =   c     + col_m) > (e =     e   + col_gep))  e = c;
	  // d is better of: score for opening a horizontal gap or score for extending a horizontal gap
	  if ((c = mCC[col] + row_m) > (d = mDD[col] + row_gep))  d = c;
	  
	  // c is score for a match
	  c = s + (this->*mMatchFunction)( row, col );
	  // put into c the best of all possible cases
	  if (e > c) c = e;
	  if (d > c) c = d;
	  
	  //--------------------------> recurrence relation <-------------------------------------------------
	  if (c <= 0)
	    {
	      c = 0;                                                  // the local alignment part
	    } else {
	      if ( c == d )                   // horizonXtal gap
		mTrace[getTraceIndex(row,col)] = TB_INSERTION;
	      else if ( c == e )              // vertical gap
		mTrace[getTraceIndex(row,col)] = TB_DELETION;
		else {                          // match
		    mTrace[getTraceIndex(row,col)] = TB_MATCH;
		}
		
	    }

	    s = mCC[col];
	    mCC[col] = c;                                              // save new score for next i
	    mDD[col] = d;
	    
	    if (mScore < c) {                                            // save maximum
	      mScore   = c;
	      mRowLast = row;
	      mColLast = col;
	    }
	}
    }
    
}   

//---------------------------------< the actual alignment algorithm >-------------------------------------------
void ImplAlignatorFullDP::performAlignmentWithoutTraceBack( const Alignandum * prow, const Alignandum * pcol) {

  register int   row, col;

  Position row_length = getRowLength();
  Position col_length = getColLength();

  Score row_gop = getRowGop();
  Score row_gep = getRowGep();
  Score col_gop = getColGop();
  Score col_gep = getColGep();

  mRowLast = 0;
  mColLast = 0;
  mScore = 0;

  Score c, e, d, s;                  // helper variables
  
  Score row_m = row_gop + row_gep;
  Score col_m = col_gop + col_gep;

  //----------------------------> Initialise affine penalty arrays <-------------------------------
    mCC[0] = 0;
    for (col = 1; col <= col_length; col++) {
      mCC[col]   = 0;
      mDD[col]   = row_gop;                               // score for horizontal gap opening
    }
 
    mCC[col_length]   = col_gop;
    
    //----------------------------> Calculate dynamic programming matrix <----------------------------
    //----------------------------> iterate over rowumns <--------------------------------------------
    for (row = 1; row <= row_length; row++) {

	// this part is different from the ordinary alignment
      s = mCC[0];
      mCC[0] = c = 0;
      
      e = col_gop;                                        // penalty for opening a vertical gap
	
      //-------------------------> iterate over cols <------------------------------------------------
      for (col = 1; col <= col_length; col++) {
 
	    //---------------------------> calculate scores <--------------------------------------------
	    // c contains score of cell above
	    // s contains score for cell [row, col-1]
	    // e is better of: score for opening a vertical gap or 
	    // score for extending a vertical gap: use col-gap-penalties
	    if ((c =   c     + col_m) > (e =     e   + col_gep))  e = c;
	    // d is better of: score for opening a horizontal gap or score for extending a horizontal gap
	    if ((c = mCC[col] + row_m) > (d = mDD[col] + row_gep))  d = c;
	    
	    // c is score for a match
	    c = s + (this->*mMatchFunction)( row, col );
	    // put into c the best of all possible cases
	    if (e > c) c = e;
	    if (d > c) c = d;
	    
      //--------------------------> recurrence relation <-------------------------------------------------
	    if (c <= 0) {
		c = 0;                                                  // the local alignment part
	    }

	    s = mCC[col];
	    mCC[col] = c;                                              // save new score for next i
	    mDD[col] = d;
	    
	    if (mScore < c) {                                            // save maximum
	      mScore   = c;
	      mRowLast = row;
	      mColLast = col;
	    }
	}
    }

}   

//------------------------------------< match functions >--------------------------------------------
Score ImplAlignatorFullDP::matchSequenceSequence( int row, int col ) {
#ifdef DEBUG
  if (verbose > LL3) {
    cout << row << " " << (int)mRowSequence[row] << "->" 
	 << col << " " << (int)mColSequence[col] << "="
	 << mSubstitutionMatrix->getScore(mRowSequence[row],mColSequence[col])
	 << endl;
  }
#endif
  return (mSubstitutionMatrix->getScore(mRowSequence[row],mColSequence[col]));
}

Score ImplAlignatorFullDP::matchSequenceSequenceTCO( int row, int col ) {
  return mSubstitutionMatrix->getScore(mRowSequence[row],mColSequence[col]);
}

Score ImplAlignatorFullDP::matchSequenceTCOSequence( int row, int col ) {
    return mSubstitutionMatrix->getScore(mRowSequence[col],mColSequence[row]);
}

Score ImplAlignatorFullDP::matchProfileSequence( int row, int col ) {
    if (mColSequence[col] == CODE_MASK)
	return MASK_VALUE;
    else
	return mRowProfile[row][mColSequence[col]];
}

Score ImplAlignatorFullDP::matchSequenceProfile( int row, int col ) {
    if (mRowSequence[row] == CODE_MASK)
	return MASK_VALUE;
    else
	return mColProfile[col][mRowSequence[row]];
}

Score ImplAlignatorFullDP::matchProfileProfile( int row, int col ) {
#ifdef DEBUG
  if (verbose > LL5) {
    cout << row << " " << col << endl;
    for (int i = 0; i < PROFILEWIDTH; i++ ) 
      cout << "("  << mRowProfile[row][i] << "/" << mColProfile[col][i] << " )" << " ";
    cout << endl;
  }
#endif

  Score score = 0;
  for (int i = 0; i < PROFILEWIDTH; i++ ) 
    score +=	mRowProfile[row][i] * mColFrequencies[col][i] + 
		mColProfile[col][i] * mRowFrequencies[row][i];
  
  return score;
}

} // namespace alignlib
