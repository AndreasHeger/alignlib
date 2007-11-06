/*
  alignlib - a library for aligning protein sequences

  $Id: ImplAlignatum.cpp,v 1.6 2004/09/16 16:02:38 aheger Exp $

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
#include <cassert>
#include "alignlib.h"
#include "AlignlibDebug.h"

#include "Alignatum.h"
#include "ImplAlignatum.h"
#include "ImplTranslator.h"
#include "Alignata.h"
#include "AlignataIterator.h"
#include "Alignandum.h"
#include "Renderer.h"
#include "HelpersTranslator.h"
#include "AlignException.h"

#ifdef WITH_DMALLOC
#include <dmalloc.h>
#endif

using namespace std;

namespace alignlib {

  //-------------------------------------------------------------------------------------------------
  /** factory functions */
  Alignatum * makeAlignatum() 
  {
    return new ImplAlignatum();
  }

  Alignatum * makeAlignatumFromString( const std::string & src, 
      Position from, 
      Position to ) 
      {
        return new ImplAlignatum( src, from, to);
      }
  Alignatum * makeAlignatumFromString( const char * src, 
      Position from, 
      Position to ) 
      {
        std::string s( src );
        return new ImplAlignatum( s, from, to);
      }

  Alignatum * makeAlignatum(const Alignandum * src, 
      const Alignata * map_this2new,
      const Position max_length) 
      {

        std::string s( src->asString() );
        Alignatum * result = new ImplAlignatum( s );

        Position length;

        if (map_this2new) {
          if (max_length == 0)
            length = map_this2new->getColTo();
          else
            length = max_length;
          result->mapOnAlignment( map_this2new, length );
        }

        return result;
      }


  //---------------------------------------------------------< constructors and destructors >--------
  ImplAlignatum::ImplAlignatum() :   
    mRepresentation(""),
    mFrom(NO_POS), 
    mTo (NO_POS),
    mLength(0), 
    mGapChar (getDefaultTranslator()->getGapChar()),
    mSeparator('\t') {
  }


  ImplAlignatum::ImplAlignatum (const std::string & representation, 
      Position from, 
      Position to) : 
        mRepresentation(representation), 
        mFrom( from), mTo( to ),
        mGapChar (getDefaultTranslator()->getGapChar()),
        mSeparator('\t') {

          mLength = mRepresentation.length();

          if (mFrom == NO_POS && mLength > 0)
            mFrom = 0; 

          if (mTo == NO_POS) 
            mTo = mFrom + mLength - countGaps();
  }

  //--------------------------------------------------------------------------------------------
  /* create a newly aligned object from ImplAlignatum, but map using ali in alignment 
   note, that the residues in the alignment are called 1..length while the residues
   in the aligned strings are 0..length-1
   */
  ImplAlignatum::ImplAlignatum (const ImplAlignatum & src ):
    mRepresentation(src.mRepresentation),
    mFrom( src.mFrom), 
    mTo ( src.mTo ),
    mLength( src.mLength), 
    mGapChar (getDefaultTranslator()->getGapChar()),
    mSeparator( src.mSeparator ) 
    {
      debug_func_cerr(5);
    }	

  //--------------------------------------------------------------------------------------------
  ImplAlignatum::~ImplAlignatum () 
    {
      debug_func_cerr(5);

    }

  //--------------------------------------------------------------------------------------------
  void ImplAlignatum::mapOnAlignment(const Alignata * map_old2new,
      const Position new_length, 
      const bool unaligned_chars) 
    {
      debug_func_cerr(5);

      std::string new_representation = "";

      // bail out on empty alignments
      if ( map_old2new->getLength() == 0)
    	  throw AlignException( "attempting to map an Alignatum object with an empty alignment ");
      
      assert( map_old2new->getRowFrom() > 0);
      // check if alignment is out-of-bounds
      assert( mLength >= map_old2new->getRowTo() );
      
      Position length = std::max( new_length, map_old2new->getColTo());

      new_representation.append( length, mGapChar );

      // get alignment start positions
      int row_from = map_old2new->getRowFrom();
      int row_to = map_old2new->getRowTo();

      // get residue numbers of terminal residues and save them in from/to
      mFrom = getResidueNumberNext(row_from); 
      mTo   = getResidueNumberPrevious(row_to-1); 

      // substitute new characters for aligned positions:
        {
          AlignataConstIterator it = map_old2new->begin();
          AlignataConstIterator it_end = map_old2new->end();

          for (; it != it_end; ++it) 
            {
              new_representation[it->mCol] = mRepresentation[it->mRow];
            }
        }

        // add unaligned characters
        if (unaligned_chars) {

          AlignataConstIterator it = map_old2new->begin();
          AlignataConstIterator it_end = map_old2new->end();

          Position last_old = it->mRow;
          Position last_new = it->mCol;

          ++it;
          // substitute new characters for aligned positions:
          for (; it != it_end; ++it) {
            Position old = it ->mRow-1;
            Position nnew = it ->mCol-1;      
            while (old - last_old > 0 && nnew - last_new > 0) {
              old--;
              nnew--;
              if (mRepresentation[old] >= 'A' &&
                  mRepresentation[old] <= 'Z' )
                new_representation[nnew] = mRepresentation[old] - 'A' + 'a';
            }
            last_old = it->mRow;
            last_new = it->mCol;
          }
        }

        mRepresentation = new_representation;
        mLength = mRepresentation.length();

    }
  //------------------------------------------------------------------------------------------------------
  ImplAlignatum * ImplAlignatum::getClone() const 
  {
    debug_func_cerr(5);

    return new ImplAlignatum(*this ); 
  }

  //-------------------------------------------------------------------------------------------------------
  ImplAlignatum * ImplAlignatum::getNew() const {
    return new ImplAlignatum();
  }

  //-------------------------------------------------------------------------------------------------------
  std::string ImplAlignatum::getString() const {
    return mRepresentation;
  }

  //-------------------------------------------------------------------------------------------------------
  const std::string & ImplAlignatum::getStringReference() const {
    return mRepresentation;
  }

  //-------------------------------------------------------------------------------------------------------
  Position ImplAlignatum::getFrom() const {
    return mFrom;
  }

  //-------------------------------------------------------------------------------------------------------
  Position ImplAlignatum::getTo() const {
    return mTo;
  }

  //-------------------------------------------------------------------------------------------------------
  void ImplAlignatum::writeRow( std::ostream & output, 
      Position segment_start, 
      Position segment_end,
      const Renderer * renderer) const {

        if (segment_start == NO_POS)
          segment_start = 0;

        if (segment_end == NO_POS || segment_end > mLength || segment_end <= segment_start) 
          segment_end = mLength;

        Position left  = getResidueNumberNext( segment_start );
        Position right = getResidueNumberPrevious( segment_end - 1);

        if (renderer) 
          output << left << getFieldSeparator() 
          << renderer->render(mRepresentation, segment_start, segment_end ) << getFieldSeparator()
          << right;
        else
          output << left << getFieldSeparator()
          << mRepresentation.substr( segment_start, segment_end - segment_start) << getFieldSeparator()
          << right;

  }

  //---------------------------------------------------------------------------------------------------
  void ImplAlignatum::readRow( std::istream & input ) {
  }

  /** get represenation */
  const std::string & ImplAlignatum::getRepresentation() const { 
    return mRepresentation; 
  } 

  /** set representation */
  void ImplAlignatum::setRepresentation( std::string & representation, 
      Position first_res,
      Position last_res) { 

        // set first residue number to 0 if not given
        if (first_res == NO_POS) 
          first_res = 0;

        mFrom = first_res;

        mRepresentation = representation; 
        mLength = mRepresentation.length();

        if (last_res == NO_POS)
          mTo = mLength - countGaps();

  }


  //---------------------------------------------------------------------------------------------------
  /** return the length of the line */
  Position ImplAlignatum::getAlignedLength() const { 
    return mLength;
  }

  //---------------------------------------------------------------------------------------------------
  /** return the length of the line without gaps */
  Position ImplAlignatum::getTrueLength() const {
    return mTo - mFrom;
  }

  //---------------------------------------------------------------------------------------------------
  /** add the specified number of gaps in the front and in the back */
  void ImplAlignatum::addGaps(int before, int after) {
    int i = 0;
    std::string x  = "";

    for (i = 0; i < before; i++) x += mGapChar;
    x+= mRepresentation;
    for (i = 0; i < after; i++)  x += mGapChar;

    mRepresentation = x;
    mLength = mRepresentation.length();
  }


  //---------------------------------------------------------------------------------------------------
  /** add one or more gaps in the middle */
  void ImplAlignatum::insertGaps( int position, Position count ) {

    std::string insertion = "";
    for (int i = 0; i < count; i++) 
      insertion += mGapChar;

    mRepresentation.insert( position, insertion );
    mLength = mRepresentation.length();
  }

  //------------------------------------------------------------------------------------------------
  /** remove leading/or trailing gaps */
  void ImplAlignatum::removeEndGaps() {
    mRepresentation.erase( 0, mRepresentation.find_first_not_of( mGapChar ));
    mRepresentation.erase( mRepresentation.find_last_not_of( mGapChar ) + 1, mRepresentation.length());
    mLength = mRepresentation.length();
  }


  //------------------------------------------------------------------------------------------------
  /** remove one or more positions from the aligned object */
  void ImplAlignatum::removeColumns( int position, Position count ) {
    mRepresentation.erase( position, position + count );
  }

  //---------------------------------------------------------------------------------------
  int ImplAlignatum::countGaps() {
    int ngaps = 0;
    Position length = mRepresentation.length();
    Position i = 0;

    for (; i < length; i++) 
      if (mRepresentation[i] == mGapChar)
        ngaps ++;

    return ngaps;
  }

  //------------------------------------------------------------------------------------------
  // Calculate the amino acid residue number of residue in string-position pos. Start counting
  // from the left
  Position ImplAlignatum::getResidueNumberNext( Position pos ) const {

    if (pos == NO_POS || pos < 0 || pos >= mRepresentation.length() ) return NO_POS;

    Position i = 0;
    // skip over terminal gaps
    while (i < pos && mRepresentation[i] == mGapChar) i++;

    Position result = mFrom;	
    for (; i < pos; i++) 
      if (mRepresentation[i] != mGapChar) 
        ++result;

    return (result);
  }
  //------------------------------------------------------------------------------------------
  // Calculate the amino acid residue number of residue in string-position pos. Start counting
  // from the right
  Position ImplAlignatum::getResidueNumberPrevious( Position pos ) const {

    if (pos == NO_POS || pos < 0 || pos >= mRepresentation.length() ) return NO_POS;

    Position i = mRepresentation.length() - 1;
    // skip over terminal gaps
    while (i >= pos && mRepresentation[i] == mGapChar) i--;

    Position result = mTo;
    for (; i > pos; --i) 
      if (mRepresentation[i] != mGapChar) 
        result--;

    return (result);
  }

  /** write into stream */
  void ImplAlignatum::write( std::ostream & output ) const {
    writeRow( output );
  }

  /** read from stream */
  void ImplAlignatum::read( std::istream & input ) {

    input >> mFrom;
    input >> mRepresentation;
    input >> mTo;

    mLength = mRepresentation.length();

  }

  //-------------------------------------------------------------------------------------------

} // namespace alignlib
