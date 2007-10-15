/*
  alignlib - a library for aligning protein sequences

  $Id: ImplRendererMView.cpp,v 1.3 2004/01/07 14:35:36 aheger Exp $

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
#include "ImplRendererMView.h"

#ifdef WITH_DMALLOC
#include <dmalloc.h>
#endif

using namespace std;

namespace alignlib {


  const TYPE_PALETTE MVIEW_COLORS[27] = {
              /*default*/ "#BEBEBE",
		  /*A*/	"#00CD00",
		  /*B*/	"#000000",
		  /*C*/	"#FF8C00",
		  /*D*/	"#2222CC",
		  /*E*/	"#2222CC",
		  /*F*/	"#228B22",
		  /*G*/	"#00CD00",
		  /*H*/	"#228B22",
		  /*I*/	"#00CD00",
		  /*J*/  "#000000",
		  /*K*/	"#CD2222",
		  /*L*/	"#00CD00",
		  /*M*/	"#00CD00",
		  /*N*/	"#A020F0",
		  /*O*/  "#000000",
		  /*P*/	"#00CD00",
		  /*Q*/	"#A020F0",
		  /*R*/	"#CD2222",
		  /*S*/	"#A020F0",
		  /*T*/	"#A020F0",
		  /*U*/	"#000000",
		  /*V*/	"#00CD00",
		  /*W*/	"#228B22",
		  /*X*/	"#474747",
		  /*Y*/	"#228B22",
		  /*Z*/	"#000000" };

/** factory functions */
Renderer * makeRendererMView( const std::string & consensus) {

  return new ImplRendererMView( &(MVIEW_COLORS[0]), consensus );
}

//---------------------------------------------------------< constructors and destructors >--------------------------------------
ImplRendererMView::ImplRendererMView ( const TYPE_PALETTE * palette, const std::string & consensus) : 
    ImplRenderer( palette ),
    mConsensus ( consensus ) {
}
		       
ImplRendererMView::~ImplRendererMView () {
}

ImplRendererMView::ImplRendererMView (const ImplRendererMView & src ) : ImplRenderer( src ), mConsensus( src.mConsensus ) {
}


//---------------------------------------------------------< methods >--------------------------------------
#define LENGTH_FACTOR 50	// number of bytes per character needed for HTML-encoding

// render the string as HTML-string with colourint according to MView-scheme
// A-Z according to scheme above, everything else is DEFAULT-COLOR
std::string ImplRendererMView::render( const std::string & representation, 			      
				       Position segment_start, 
				       Position segment_end ) const {
    

  long estimated_size = LENGTH_FACTOR * representation.size();

  debug_cerr( 5, "Trying to allocate " << estimated_size << " bytes of memory for ImplRendererMView::render()" );

  char * buffer = new char[estimated_size];

  int lastcolor = 0;
  int color_code = 0;
  
  sprintf( buffer, "<FONT COLOR=\"%s\">", mPalette[color_code]);
  long pos = 22;
 
  for (Position i = segment_start; i < segment_end; i++) {
    
    char c = representation[i];

    if (c > 'Z' || c < 'A' || c != mConsensus[i]) 
      color_code = 0;
    else 
      color_code = c - 'A' + 1;
    
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

