/*
  alignlib - a library for aligning protein sequences

  $Id: Renderer.h,v 1.3 2004/03/19 18:23:41 aheger Exp $

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

#ifndef RENDERER_H
#define RENDERER_H 1

#include <string>
#include "alignlib.h"

namespace alignlib {

typedef char TYPE_PALETTE[10];
 
/** @short Interface definition of Renderer objects.
   
    Renderers are a plugin for multiple alignments. It takes a string representation
    of a row in the multiple alignment sequence together with the coordinates and returns 
    a new string.
    
    A common usage is to color the multiple alignment by inserting HTML statements 
    between the residues.
    
    This class is a protocol class and as such defines only the general interface.
    
    @author Andreas Heger
    @version $Id: Renderer.h,v 1.3 2004/03/19 18:23:41 aheger Exp $
*/
class Renderer {

 public:
    // constructors and desctructors

    /** default constructor */
    Renderer  ();
    
    /** copy constructor */
    Renderer  (const Renderer &);

    /** destructor */
    virtual ~Renderer ();

    // render one string into a different representation
    virtual std::string render( const std::string & representation, 			      
				Position segment_start, 
				Position segment_end ) const = 0;


};

 
}

#endif /* RENDERER_H */

