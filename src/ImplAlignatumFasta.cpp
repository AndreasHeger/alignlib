/*
  alignlib - a library for aligning protein sequences

  $Id: ImplAlignatumFasta.cpp,v 1.2 2004/01/07 14:35:34 aheger Exp $

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
#include "alignlib_fwd.h"
#include "alignlib_interfaces.h"
#include "AlignlibDebug.h"
#include "ImplAlignatumFasta.h"
#include "ImplTranslator.h"
#include "Alignment.h"
#include "AlignmentIterator.h"
#include "Alignandum.h"
#include "Renderer.h"
#include "HelpersTranslator.h"

#ifdef WITH_DMALLOC
#include <dmalloc.h>
#endif

using namespace std;

namespace alignlib 
{

	//--------------------------------------------------------------------------------------------
  	/** factory functions */
  	HAlignatum makeAlignatumFasta() 
  	{
  		return HAlignatum( new ImplAlignatumFasta() );
  	}

  	HAlignatum makeAlignatumFasta( const std::string & description,
				  const std::string & src, 
				  Position from, 
				  Position to ) 
  	{
  		return HAlignatum( new ImplAlignatumFasta( description, src, from, to) );
  	}


//---------------------------------------------------------< constructors and destructors >--------
ImplAlignatumFasta::ImplAlignatumFasta() :   
  ImplAlignatum(),
  mDescription("") {
}

ImplAlignatumFasta::ImplAlignatumFasta (const std::string & description,
					const std::string & representation, 
					Position from, 
					Position to ) :
					
  ImplAlignatum( representation, from, to ),
  mDescription( description ) {
}

//--------------------------------------------------------------------------------------------
/* create a newly aligned object from ImplAlignatumFasta, but map using ali in alignment 
   note, that the residues in the alignment are called 1..length while the residues
   in the aligned strings are 0..length-1
*/
ImplAlignatumFasta::ImplAlignatumFasta (const ImplAlignatumFasta & src ):
  ImplAlignatum(src), 
  mDescription(src.mDescription) 
{
  debug_func_cerr(5);

}	
	       
//--------------------------------------------------------------------------------------------
ImplAlignatumFasta::~ImplAlignatumFasta () 
{
  debug_func_cerr(5);
}

//------------------------------------------------------------------------------------------------------
HAlignatum ImplAlignatumFasta::getClone() const 
{
	debug_func_cerr(5);
	return HAlignatum( new ImplAlignatumFasta(*this ) ); 
}

//-------------------------------------------------------------------------------------------------------
HAlignatum ImplAlignatumFasta::getNew() const 
{
	debug_func_cerr(5);
	return HAlignatum( new ImplAlignatumFasta() );
}


//-------------------------------------------------------------------------------------------------------
void ImplAlignatumFasta::write( std::ostream & output,
		const HRenderer & renderer,
		Position segment_start, 
		Position segment_end ) const 
		{
  
  output << setw(30) << mDescription;

  ImplAlignatum::write( output, renderer, segment_start, segment_end );
}

/** write into stream */
void ImplAlignatumFasta::write( std::ostream & output ) const 
{
  output << ">" << mDescription << endl << getRepresentation() << endl;
}

/** read from stream */
void ImplAlignatumFasta::read( std::istream & input ) 
{

#define MAX_CHUNK 10000

  char * buffer = new char[MAX_CHUNK];
  
  while ( (input.peek() != '>') && 
	  !input.eof() ) {
    input.getline( buffer, MAX_CHUNK);
  }

  if (input.eof())
    return;
    
  input.get();
  input.getline(buffer, MAX_CHUNK);
    
  mDescription = buffer;

  // erase end of line
  mDescription.erase(mDescription.size(),1); 
  
  std::string representation("");
   
  // build the sequence character-wise
  while ( (input.peek() != '>') && 
	  !input.eof() ) {
      
    input.getline( buffer, MAX_CHUNK);
      
    for (unsigned int i = 0; i < strlen(buffer); i++) 
      if (getDefaultTranslator()->isValidChar( buffer[i] )) 
	representation += buffer[i];
  }
  
  setRepresentation(representation);

  delete [] buffer;

}

//-------------------------------------------------------------------------------------------

} // namespace alignlib






