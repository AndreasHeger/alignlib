/*
  alignlib - a library for aligning protein sequences

  $Id: ImplRendererColumn.h,v 1.2 2004/01/07 14:35:36 aheger Exp $

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


#if HAVE_CONFIG_H
#include <config.h>
#endif

#ifndef IMPL_RENDERER_COLUMN_H
#define IMPL_RENDERER_COLUMN_H 1

#include "alignlib_fwd.h"
#include "ImplRendererPalette.h"

namespace alignlib 
{

 /** render one string into another
     0 = default color
   
      @author Andreas Heger
      @version $Id: ImplRendererColumn.h,v 1.2 2004/01/07 14:35:36 aheger Exp $
      @short protocol class for log-odders
      
  */
    
class ImplRendererColumn : public ImplRendererPalette 
{

 public:
    // constructors and desctructors

    /** default constructor */
    ImplRendererColumn( const TYPE_PALETTE * palette, const unsigned char * color_code );
    
    /** copy constructor */
    ImplRendererColumn  (const ImplRendererColumn &);

    /** destructor */
    virtual ~ImplRendererColumn ();

    // render one string into a different representation
    virtual std::string render( const std::string & representation, 			      
				Position segment_start, 
				Position segment_end ) const ;
    
 private:
    /** the consensus string according which to color */
    const unsigned char * mColorCode;


};

 
}

#endif /* IMPL_RENDERER_MVIEW_H */

