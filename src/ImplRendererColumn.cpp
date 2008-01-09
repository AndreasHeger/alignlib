/*
  alignlib - a library for aligning protein sequences

  $Id: ImplRendererColumn.cpp,v 1.2 2004/01/07 14:35:36 aheger Exp $

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
#include <string>
#include <stdio.h>
#include "alignlib.h"
#include "AlignlibDebug.h"
#include "Renderer.h"
#include "ImplRendererColumn.h"
#include "HelpersRenderer.h"


using namespace std;

namespace alignlib {

/** factory functions */
HRenderer makeRendererColumn( const unsigned char * color_code, const TYPE_PALETTE * palette ) 
{
    if (!palette) 
    	palette = getDefaultPalette();
	    
    return HRenderer( new ImplRendererColumn( palette, color_code ) );
}
	
//---------------------------------------------------------< constructors and destructors >--------------------------------------
ImplRendererColumn::ImplRendererColumn ( const TYPE_PALETTE * palette, const unsigned char * color_code ) : 
    ImplRendererPalette( palette ),
    mColorCode( color_code ) {
}
		       
ImplRendererColumn::~ImplRendererColumn () {
}

ImplRendererColumn::ImplRendererColumn (const ImplRendererColumn & src ) : 
    ImplRendererPalette( src ), 
    mColorCode( src.mColorCode ) 
    {
}

#define LENGTH_FACTOR 50	// number of bytes per character needed for HTML-encoding

// render the string as HTML-string with colourint according to Column-scheme
// A-Z according to scheme above, everything else is DEFAULT-COLOR
std::string ImplRendererColumn::render( const std::string & representation, 			      
					Position segment_start, 
					Position segment_end ) const {
    
  long estimated_size = LENGTH_FACTOR * representation.size();

  debug_cerr( 5, "Trying to allocate " << estimated_size << " bytes of memory for ImplRendererColumn::render()" );
  
  char * buffer = new char[estimated_size];
  
  int lastcolor = 0;
  unsigned char color_code = 0;

  sprintf( buffer, "<FONT COLOR=\"%s\">", mPalette[color_code]);
  long pos = 22;

  for (Position i = segment_start; i < segment_end; i++) {
    
    char c = representation[i];

    if (c > 'Z' || c < 'A') 
      color_code = 0;
    else 
      color_code = mColorCode[i];
    
    if (lastcolor != color_code) {
      sprintf( &buffer[pos], "</FONT><FONT COLOR=\"%s\">",  mPalette[color_code]);
      lastcolor = color_code;
      pos += 29;
    }
    
    buffer[pos++] = c;
  }
  
  sprintf(&buffer[pos], "</FONT>");
  buffer[pos+7] = '\0';
  
  std::string result( buffer );
  delete [] buffer;
    
  return result;
}

} // namespace alignlib
