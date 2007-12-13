/*
  alignlib - a library for aligning protein sequences

  $Id: ImplAlignatumMView.cpp,v 1.3 2004/03/19 18:23:41 aheger Exp $

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
#include "ImplAlignatum.h"
#include "ImplTranslator.h"
#include "Alignata.h"
#include "AlignataIterator.h"
#include "Alignandum.h"

#ifdef WITH_DMALLOC
#include <dmalloc.h>
#endif

using namespace std;

namespace alignlib {

    /** factory functions */
  Alignatum * makeAlignatum(const Alignatum * src, 
			    const Alignata * ali, 
			    bool skip_gaps,
			    bool is_in_row) {
    return new ImplAlignatum( src, ali, skip_gaps, is_in_row);
  }

  extern const alignlib::ImplTranslator DEFAULT_TRANSLATOR; // defined in ImplTranslatorBlosum.cpp

//---------------------------------------------------------< constructors and destructors >--------------------------------------
ImplAlignatum::ImplAlignatum() :   
  mRepresentation(""),
  mFrom(0), 
  mTo (0),
  mLength(0), 
  mGapChar (DEFAULT_TRANSLATOR.getGapChar()) {
}


ImplAlignatum::ImplAlignatum (const std::string & representation, 
			      Position from, 
			      Position to) : 
  mRepresentation(representation), 
  mFrom( from), mTo( to ),
  mGapChar (DEFAULT_TRANSLATOR.getGapChar()) {

  mLength = mRepresentation.length();
  
  if (mTo == 0) 
    mTo = mFrom + mLength - countGaps() - 1;
}

//--------------------------------------------------------------------------------------------
/* create a newly aligned object from ImplAlignatum, but map using ali in alignment 
   note, that the residues in the alignment are called 1..length while the residues
   in the aligned strings are 0..length-1
*/
ImplAlignatum::ImplAlignatum (const ImplAlignatum & src, 
			      const Alignata * ali, 
			      bool skip_gaps, 
			      bool is_in_row ) : 
  mRepresentation(src.mRepresentation),
  mFrom( src.mFrom), 
  mTo ( src.mTo ),
  mLength( src.mLength), 
  mGapChar (DEFAULT_TRANSLATOR.getGapChar()) 
{
  debug_func_cerr(5);

    
    if (ali) 
	mapOnAlignment( ali, skip_gaps, is_in_row);
}	
	       
//--------------------------------------------------------------------------------------------
ImplAlignatum::~ImplAlignatum () 
{
  debug_func_cerr(5);

}

//--------------------------------------------------------------------------------------------
ImplAlignatum::ImplAlignatum (const ImplAlignatum & src ) : 
    mRepresentation(src.mRepresentation),
    mFrom (src.mFrom) , 
    mTo (src.mTo), 
    mLength( src.mLength ),
    mGapChar (src.mGapChar) {
}

//--------------------------------------------------------------------------------------------
ImplAlignatum::ImplAlignatum(const Alignandum * src, 
			     const Alignata * ali, 
			     bool skip_gaps,
			     bool is_in_row) : 

    mRepresentation(src->asString()), 
    mFrom(1), 
    mTo(0),
    mLength(),
    mGapChar (DEFAULT_TRANSLATOR.getGapChar()) {
    
    mLength = mRepresentation.length();
    mTo = mFrom + mLength - countGaps() - 1;
    
    if (ali) 
	mapOnAlignment( ali, skip_gaps, is_in_row);

}
    

  
//--------------------------------------------------------------------------------------------
void ImplAlignatum::mapOnAlignment(const Alignata * ali, bool skip_gaps, bool is_in_row ) 
{
  debug_func_cerr(5);


  std::string new_representation = "";

  // get alignment start positions
  int row_from = ali->getRowFrom();
  int col_from = ali->getColFrom();

  int row_to   = ali->getRowTo();
  int col_to   = ali->getColTo();

  // get amino-acid start positions and save them in from/to
  if (is_in_row) {
      mFrom = getResidueNumbernext(row_from-1), mTo = getResidueNumberprevious(row_to-1); 
  } else {
      mFrom = getResidueNumbernext(col_from-1), mTo = getResidueNumberprevious(col_to-1); 
  }

  int d;
  int last_row = row_from;
  int last_col = col_from;
  int current_row, current_col;
  
  AlignataConstIterator it = ali->begin();
  AlignataConstIterator it_end = ali->end();

#ifdef DEBUG
  if (is_in_row) 
    cout << "Is in row" << endl;
  else 
    cout << "Is in col" << endl;

  { 
    for (unsigned int i = 0; i < mRepresentation.length(); i++) 
      cout << i << " " << mRepresentation[i] << endl;
  }
#endif

  //--> iterate over aligned pairs 
  for (; it != it_end; ++it) {
    ResiduePAIR pair = *it;
    
    // if no difference in either sequence
    current_row = pair.mRow;
    current_col = pair.mCol;

    if (is_in_row) {
      if (!skip_gaps) {						// if no collapse add insertions
	while (++last_row < current_row) 
	  new_representation += mRepresentation[last_row-1];
      }
      
      if ( (d = (current_col - last_col)) > 1 ) 		// add deletions
	while (--d > 0)
	  new_representation += mGapChar;
      new_representation += mRepresentation[current_row-1];			// add match
    } else {
      if (!skip_gaps) {						// if no collapse add insertions
	while (++last_col < current_col) 
	  new_representation += mRepresentation[last_col-1];
      }
      
      if ( (d = (current_row - last_row)) > 1 ) 		// add deletions
	while (--d > 0)
	  new_representation += mGapChar;
      
      new_representation += mRepresentation[current_col-1];			// add match
    }
    
    last_row = current_row;
    last_col = current_col;
  }
  
  mRepresentation = new_representation;
  mLength = mRepresentation.length();
  
}
//--------------------------------------------------------------------------------------------------------------------------------
ImplAlignatum * ImplAlignatum::getClone( const Alignata * ali, bool skip_gaps, bool is_in_row) const 
{
  debug_func_cerr(5);

  return new ImplAlignatum(*this, ali, skip_gaps, is_in_row ); 
}

//--------------------------------------------------------------------------------------------------------------------------------
ImplAlignatum * ImplAlignatum::getNew() const {
  return new ImplAlignatum();
}

//--------------------------------------------------------------------------------------------------------------------------------
std::string ImplAlignatum::getString() const {
    return mRepresentation;
}

//--------------------------------------------------------------------------------------------------------------------------------
const std::string & ImplAlignatum::getStringReference() const {
    return mRepresentation;
}
       
//--------------------------------------------------------------------------------------------------------------------------------
Position ImplAlignatum::getFrom() const {
    return mFrom;
}
    
//--------------------------------------------------------------------------------------------------------------------------------
Position ImplAlignatum::getTo() const {
    return mTo;
}

//--------------------------------------------------------------------------------------------------------------------------------
void ImplAlignatum::writeRow( std::ostream & output, 
			      Position segment_start, 
			      Position segment_end) const {

    if (segment_end >= mLength || segment_end <= segment_start) 
	segment_end = mLength;
  
    Position left  = getResidueNumbernext( segment_start );
    Position right = getResidueNumberprevious( segment_end - 1);

    output << setw(6) << left << " " 
	   << mRepresentation.substr( segment_start, segment_end - segment_start) 
	   << setw(6) << right;
}

//--------------------------------------------------------------------------------------------------------------------------------
void ImplAlignatum::readRow( std::istream & input ) {
}

//--------------------------------------------------------------------------------------------------------------------------------
/** return the length of the line */
Position ImplAlignatum::getAlignedLength() const { 
    return mLength;
}

//--------------------------------------------------------------------------------------------------------------------------------
/** return the length of the line without gaps */
Position ImplAlignatum::getFullLength() const {
    return mTo - mFrom + 1;
}
    
//--------------------------------------------------------------------------------------------------------------------------------
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

//--------------------------------------------------------------------------------------------------------------------------------
/** add one or more gaps in the middle */
void ImplAlignatum::insertGaps( int position, Position count ) {

  std::string insertion = "";
  for (int i = 0; i < count; i++) 
    insertion += mGapChar;

  mRepresentation.insert( position, insertion );
  mLength = mRepresentation.length();
}

//--------------------------------------------------------------------------------------------------------------------------------
/** remove leading/or trailing gaps */
void ImplAlignatum::removeEndGaps() {
  mRepresentation.erase( 0, mRepresentation.find_first_not_of( mGapChar ));
  mRepresentation.erase( mRepresentation.find_last_not_of( mGapChar ) + 1, mRepresentation.length());
  mLength = mRepresentation.length();
}


//--------------------------------------------------------------------------------------------------------------------------------
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
Position ImplAlignatum::getResidueNumbernext( Position pos ) const {

    Position i = 0;
    Position result = mFrom;
    while (i < pos && mRepresentation[i] == mGapChar) i++;
    
    for (; i < pos; i++) 
	if (mRepresentation[i] != mGapChar) 
	    result++;
    
  return (result);
}
//------------------------------------------------------------------------------------------
// Calculate the amino acid residue number of residue in string-position pos. Start counting
// from the right
Position ImplAlignatum::getResidueNumberprevious( Position pos ) const {
// #ifdef DEBUG
//     if (verbose > LL3) cout << "Position ImplAlignatum::getResidueNumberprevious( Position) const " << endl;
// #endif

    Position i = mRepresentation.length() - 1;
    Position result = mTo;
    while (i >= pos && mRepresentation[i] == mGapChar) i--;
    
    for (; i > pos; i--) 
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

};


//--------------------------------------------------------------------------------------------------------------------------------

} // namespace alignlib
