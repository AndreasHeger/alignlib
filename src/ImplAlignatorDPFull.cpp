/*
  alignlib - a library for aligning protein sequences

  $Id: ImplAlignatorDPFull.cpp,v 1.1 2005/02/24 11:08:24 aheger Exp $

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

#include <sys/types.h>
#include <sys/socket.h>
#include <sys/time.h>      


#include "SubstitutionMatrix.h"
#include "HelpersSubstitutionMatrix.h"
#include "ImplSubstitutionMatrixAA.h"

#include "Alignata.h"
#include "HelpersAlignata.h"

#include "Alignandum.h"
#include "ImplAlignatorDPFull.h"
#include "Alignator.h"

//implementation of Alignandum-objects
#include "ImplSequence.h"
#include "ImplProfile.h"

#include "ImplTranslator.h" // for direct access to mask_code


#ifdef WITH_DMALLOC
#include <dmalloc.h>
#endif

using namespace std;

namespace alignlib
{

  Alignator * makeAlignatorDPFull( AlignmentType alignment_type,
      Score gop, Score gep,
      bool penalize_left,
      bool penalize_right )
      {
        return new ImplAlignatorDPFull( alignment_type, gop, gep, gop, gep,
            penalize_left, penalize_right, penalize_left, penalize_right);
      }

  /* How to write a fast algorithm:
     My design objective here was to not duplicate the algorithmic code without penalizing too much
     for function indirection. It the way I do it below, there is one indirection for every call to
     calculate a match function, except for sequence-sequence comparisons, where there would have been
     two.

   */

  //<--------------< Global functions and pointers for fast determination of match score. >-------
  // I don't know how to use member functions as function pointers. After all, this is what inheritence is for. A possibility would be
  // to automatically subclass an alignator-object.However, I do not like this idea, since this assumes that the parent has information 
  // about the child. On the other hand, via the inlining mechanism a function call could be saved. The problem is when you want to change
  // the algorithm by overloading align. Then the parent functions () do not know, what the child is. Therefore I use the static 
  // functions. Maybe it is possible to separate the algorithm and the type-decision into different classes that interact.
  //
  // The danger of static functions is that the global pointers are unsafe, i.e. there exist just one copy for all alignator-objects, and
  // this code will never be threadsafe.
  //----------------------------------------------------------------------------------------------

  ImplAlignatorDPFull::ImplAlignatorDPFull( AlignmentType alignment_type, 
      Score row_gop, Score row_gep, 
      Score col_gop, Score col_gep,
      bool penalize_row_left, bool penalize_row_right,
      bool penalize_col_left, bool penalize_col_right) :
        ImplAlignatorDP( alignment_type, row_gop, row_gep, col_gop, col_gep,
            penalize_row_left, penalize_row_right, penalize_col_left, penalize_col_right ), 
            mTraceMatrix(NULL), mTraceRowStarts(NULL), 
            mRowFrom(NO_POS), mRowTo( NO_POS),
            mRowLast(NO_POS), mColLast(NO_POS)
            {
            }

  //----------------------------------------------------------------------------------------------
  ImplAlignatorDPFull::ImplAlignatorDPFull( const ImplAlignatorDPFull & src ) : ImplAlignatorDP( src ) 
  {
    debug_func_cerr(5);


    mTraceMatrix = NULL;
    mTraceRowStarts = NULL;
    mRowFrom = NO_POS;
    mRowTo = NO_POS;
    mRowLast = NO_POS;
    mColLast = NO_POS;

    mPenalizeRowLeft = src.mPenalizeRowLeft;
    mPenalizeRowRight = src.mPenalizeRowRight;
    mPenalizeColLeft = src.mPenalizeColLeft;
    mPenalizeColRight = src.mPenalizeColRight;

  }

  //------------------------------------------------------------------------------------------------
  ImplAlignatorDPFull::~ImplAlignatorDPFull() 
    {
      debug_func_cerr(5);


    }

  //------------------------------------------------------------------------------------------------
  void ImplAlignatorDPFull::startUp(const Alignandum * row, const Alignandum *col, Alignata * ali)

    {
      debug_func_cerr(5);

      ImplAlignatorDP::startUp(row, col, ali);  

      mRowLast = NO_POS;
      mColLast = NO_POS;

      //---------------------------------------------
      // setup Traceback-matrix
      // - matrix can have variable length per row
      // - there is a 0 element in each col (col_front - 1). 
      // - in other words: indexing for columns and rows starts at 0 in the trace matrix, but there is
      //   a -1 element for each column
      debug_cerr( 5, "allocating start positions for " << mIterator->row_size() << " rows." );
      mTraceRowStarts = new TraceIndex[mIterator->row_size() + 1];
      // shift index one up, so that there is a -1 element    
      ++mTraceRowStarts;
      mTraceRowStarts[-1] = 0;
      TraceIndex matrix_size =  1 + mIterator->col_size();

      Iterator2D::const_iterator rit(mIterator->row_begin()), rend(mIterator->row_end());
      mRowFrom = *rit;
      mRowTo = *rend;

      for (unsigned row = 0; rit != rend; ++rit, ++row)
        {
          mTraceRowStarts[row] = matrix_size;
          matrix_size += 1 + mIterator->col_size( *rit );    	
          debug_cerr( 6, "trace matrix row=" << row << " matrix_size=" << matrix_size << " col_size=" << mIterator->col_size(*rit));
        }

      mMatrixSize = matrix_size;
      debug_cerr( 5, "allocating trace matrix for " << matrix_size << " elements." );
      mTraceMatrix = new TraceEntry[matrix_size];
      for (TraceIndex i = 0; i < matrix_size; i++)
        mTraceMatrix[i] = TB_STOP;

    }

  //--------------------------------------------------------------------------------------------------------------
  void ImplAlignatorDPFull::cleanUp(const Alignandum * row, const Alignandum *col, Alignata * ali)
    {
      if (mTraceMatrix != NULL) { delete [] mTraceMatrix; mTraceMatrix = NULL; }
      if (mTraceRowStarts != NULL) { --mTraceRowStarts; delete [] mTraceRowStarts; mTraceRowStarts = NULL; }    

      ImplAlignatorDP::cleanUp(row, col, ali);
    }       

  //-------------------------------------< BackTrace >----------------------------------------------------------------------

  // wrapping around for col but not for row, because otherwise there could be an infinite loop.
#define PREVCOL { if (--col < 0) col = mColLength - 1; }  

  void ImplAlignatorDPFull::traceBack( const Alignandum * prow, const Alignandum * pcol, Alignata * result) 
    {	
      debug_func_cerr(5);

      int t;

      if (mRowLast == NO_POS || mColLast == NO_POS)
        return;

#ifdef DEBUG
        {
          debug_cerr(5, "Trace matrix" );

            {
            	debug_cerr_start( 5, setw(6) << " " );
            	for (Position c = -1; c < mIterator->col_size(); ++c) 
            		debug_cerr_add( 5, setw(4) << c );
            	debug_cerr_add( 5, std::endl );
            }

            Iterator2D::const_iterator rit(mIterator->row_begin()), rend(mIterator->row_end());
            for (; rit != rend; ++rit)
              {
                Position row = *rit;
                debug_cerr_start( 5, setw(6) << row );
                for (Position col = 0; col < mIterator->col_front(row) - 1; ++col)
                	debug_cerr_add( 5, setw(4) << "-" );
                for (Position col = mIterator->col_front(row) - 1; col <= mIterator->col_back(row); ++col)
                  {
                	debug_cerr_add( 5, setw(4) );
                    switch (mTraceMatrix[getTraceIndex(row,col)])
                    {
                    case TB_STOP:       debug_cerr_add( 5, "o" ); break;
                    case TB_DELETION:   debug_cerr_add( 5, "<" ); break;
                    case TB_INSERTION:  debug_cerr_add( 5, "^" ); break;
                    case TB_MATCH:      debug_cerr_add( 5, "=" ); break;
                    case TB_WRAP:       debug_cerr_add( 5, "@" ); break;
                    default: debug_cerr_add( 5, "#" ); break;
                    }
                  }
                for (Position col = mIterator->col_back(row) + 1; col < mIterator->col_size(); ++col)
                  debug_cerr_add( 5, setw(4) << "-" );                  
                debug_cerr_add( 5, std::endl );
              }
            debug_cerr( 5, "traceback starts in cell ("<< mRowLast << "," << mColLast << ") with score " << mScore );
        }

#endif

        Position col = mColLast;
        Position row = mRowLast;
        int ngaps = 0;
        
        Position row_from = mIterator->row_front();
        
        t = mTraceMatrix[getTraceIndex(row,col)];

        while ( t != TB_STOP ) 
          {

            switch (t) {
            case TB_DELETION :
              col--;
              ngaps++;
              if (col < 1) 
                row--;
              break;
            case TB_INSERTION :
              ngaps++;
              row--;
              break;
            case TB_MATCH :
              result->addPair( new ResiduePAIR( row, col, mScorer->getScore( row, col)));
              row--;
              col--;
              break;
            case TB_WRAP :
              col = mIterator->col_back( row );
              break;
            default:
              throw AlignException("Unknown matrix command in TraceBack");
              break;
            }
            if (row < row_from) break;
            t = mTraceMatrix[getTraceIndex(row,col)];
          }
        result->setScore ( mScore );
    }  

  //---------------------------------< the actual alignment algorithm >-------------------------------------------

  //-----------------------------------------------------------------------------------  
  void ImplAlignatorDPFull::performAlignment( const Alignandum * prow, const Alignandum * pcol, Alignata * ali)
    {

      debug_func_cerr(5);
      
      Score row_gop = getRowGop();
      Score row_gep = getRowGep();
      Score col_gop = getColGop();
      Score col_gep = getColGep();

      Score c, e, d, s;                  // helper variables
      c = e = d = s = 0;

      Score row_m = row_gop + row_gep;
      Score col_m = col_gop + col_gep;

      /*
      ------> col
      | CC-> 
      | DD->
      |
      |
      row

      For each cell:
        s   mCC/mDD   
         \  |
          \ |
      c/e-- x

      c/mCC: last op was match
      e/mDD: last op was gap
       */

      //---> Initialise affine penalty arrays <--------------------------

      //-----------------------------------------------------------------
        {

          Position row = mIterator->row_front();
          Iterator2D::const_iterator cit(mIterator->col_begin(row)), cend(mIterator->col_end(row));
          assert( (*cit) -1 >= -1);
          mCC[(*cit)-1] = 0;

          switch (mAlignmentType)
          {
          case ALIGNMENT_LOCAL:
            for (; cit != cend; ++cit)
              {
                Position col = *cit;	    
                mCC[col]   = 0;
                mDD[col]   = row_gop;                               // score for horizontal gap opening
              }
            mCC[mIterator->col_back()] = col_gop;
            break;
          case ALIGNMENT_WRAP:      
            for (; cit != cend; ++cit)
              {
                Position col = *cit;	    
                mCC[col]   = 0;
                mDD[col]   = row_gop;                               // score for horizontal gap opening
              }
            mCC[mIterator->col_back()] = col_gop;
            break;
          case ALIGNMENT_GLOBAL:
            /* set initial values for upper border */      
            if (mPenalizeRowLeft)
              {
                for (; cit != cend; ++cit)
                  {
                    Position col = *cit;	    
                    mCC[col]   = row_gop + row_gep * (col+1);
                    mDD[col]   = mCC[col];                              // add score for gap opening
                  }
              }
            else
              {
                for (; cit != cend; ++cit)
                  {
                    Position col = *cit;	    		
                    mCC[col]   = 0;
                    mDD[col]   = row_gop;                               // add score for gap opening
                  }
              }
          }
        }

        //----------------------------> Calculate dynamic programming matrix <----------------------------
        //----------------------------> iterate over rowumns <--------------------------------------------
        Iterator2D::const_iterator rit(mIterator->row_begin()), rend(mIterator->row_end());

        debug_cerr( 5, "aligning within the following coordinates: row= "
            << *mIterator->row_begin() << "-" <<  *mIterator->row_end() << ":" << mIterator->row_size() << " col=" 
            << *mIterator->col_begin() << "-" <<  *mIterator->col_end() << ":" << mIterator->col_size() );
        debug_cerr( 5, "penalties=" << mPenalizeRowLeft << " " << mPenalizeRowRight 
            << " " << mPenalizeColLeft << " " << mPenalizeColRight );

        for (; rit != rend; ++rit)
          {
            Position row = *rit;
            Position row_length = mIterator->row_size();
            Position row_from   = mIterator->row_front( row );
            Position col_length = mIterator->col_size( row );

            Iterator2D::const_iterator cit(mIterator->col_begin(row)), cend(mIterator->col_end(row));
            Position col_from = *cit;

            switch (mAlignmentType)
            {
            case ALIGNMENT_LOCAL:
              s = mCC[col_from - 1];
              mCC[col_from - 1] = c = 0;
              e = col_gop;  // penalty for opening a vertical gap
              break;
            case ALIGNMENT_WRAP:
              // the wrapping around part	  
              if (mCC[col_length-1] > 0)
                {
                  mCC[col_from - 1] = c = mCC[col_length-1];
                  mTraceMatrix[getTraceIndex(row-1,-1)] = TB_WRAP;
                }
              else
                {
                  mCC[col_from - 1] = c = 0;
                }

              s = mCC[col_from - 1];
              e = col_gop; // penalty for opening a vertical gap
              break;
            case ALIGNMENT_GLOBAL:
              /* set initial values for left border */
              // does this work with iterator2D?
              if (mPenalizeColLeft)
                {
                  s = mCC[col_from-1];         
                  mCC[col_from-1] += col_gep;
                  if (row == row_from)
                    mCC[col_from-1] += col_gop;
                  e = c = col_gop + col_gep * (row + 1);                  
                }
              else
                {
                  s = 0;
                  c = 0;
                  e = col_gop;
                }
              break;
            }

            //-------------------------> iterate over cols <------------------------------------------------	
            for (; cit != cend; ++cit)
              {
                Position col = *cit;

                //---------------------------> calculate scores <--------------------------------------------
                // c contains score of cell left
                // s contains score for cell [row-1, col-1]
                // e is better of: score for opening a horizontal gap or score for extending a horizontal gap: 
                // use col-gap-penalties
                bool new_horizontal_gap = false;
                if ((c = c + col_m) > (e = e + col_gep))  
                  {
                    new_horizontal_gap = true;
                    e = c;
                  }       
                // d is better of: score for opening a vertical gap or score for extending a vertical gap
                bool new_vertical_gap = false;
                if ((c = mCC[col] + row_m) > (d = mDD[col] + row_gep))
                  {
                    new_vertical_gap = true;
                    d = c;
                  }

                // c is score for a match
                c = s + mScorer->getScore(row,col);

                // put into c the best of all possible cases
                if (e > c) c = e;
                if (d > c) c = d;

                //--------------------------> recurrence relation <-------------------------------------------------
                if (mAlignmentType == ALIGNMENT_LOCAL && c <= 0)
                  {
                    c = 0;                                                  // the local alignment part
                  }
                  else
                  {
                    if ( c == d )                   // vertical gap
                      {
                        mTraceMatrix[getTraceIndex(row,col)] = TB_INSERTION;
                        // make sure that if we backtrack via gap extension
                        // if this is an old gap.
                        if (!new_vertical_gap)
                           mTraceMatrix[getTraceIndex(row-1,col)] = TB_INSERTION;
                      }
                    else if ( c == e )              // horizontal gap
                      {
                        mTraceMatrix[getTraceIndex(row,col)] = TB_DELETION;
                        if (!new_horizontal_gap)
                          mTraceMatrix[getTraceIndex(row,col-1)] = TB_DELETION;
                      }
                    else                           // match
                      mTraceMatrix[getTraceIndex(row,col)] = TB_MATCH;
                  }

                debug_cerr( 5, " row=" << row << " col=" << col << " c=" << c << " e=" << e << " d=" << d << " s=" << s
                        << " mCC=" << mCC[col] << " mDD=" << mDD[col]
                        << " index=" << getTraceIndex(row,col) <<  " mScore=" << mScore << " : "                               
                        << (char*) (( c == d ) ? "insertion" : (( c == e ) ? "deletion" : "match") ) );      

                s = mCC[col];
                mCC[col] = c;                                              // save new score for next i
                mDD[col] = d;

                /* retrieve maximum score. If penalty has to be paid for right end-gaps, 
	       then restrict search correspondingly;
                 */
                if (mAlignmentType == ALIGNMENT_GLOBAL)
                  {                    
                    if (mPenalizeRowRight && row < row_length - 1) 
                      continue;
                    if (mPenalizeColRight && col < col_length - 1) 
                      continue;
                  }

                if (mScore < c)
                  {
                    // save maximum
                	debug_cerr( 5, " updating score from cell " << row << "," << col << " with " << c );
                    mScore   = c;
                    mRowLast = row;
                    mColLast = col;
                  }
              }
          }

        traceBack(prow, pcol, ali);
    }

} // namespace alignlib
