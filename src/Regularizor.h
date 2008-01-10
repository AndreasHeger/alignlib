/*
  alignlib - a library for aligning protein sequences

  $Id: Regularizor.h,v 1.2 2004/01/07 14:35:37 aheger Exp $

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

#ifndef REGULARIZOR_H
#define REGULARIZOR_H 1

#include "alignlib.h"

namespace alignlib 
{
	
  /** Protocoll class for objects, that regularize profile columns.
      
      This class is a protocol class and as such defines only the general interface.
      
      @author Andreas Heger
      @version $Id: Regularizor.h,v 1.2 2004/01/07 14:35:37 aheger Exp $
      @short protocol class for sequence weighters
      
  */

class Regularizor 
{
 public:
    // constructors and desctructors

    /** default constructor */
    Regularizor  ();
    
    /** copy constructor */
    Regularizor  (const Regularizor &);

    /** destructor */
    virtual ~Regularizor ();
    
    /** copy the counts into the frequencies and regularize them by doing so. */
    virtual void fillFrequencies( 
    		FrequencyMatrix * frequencies, 
    		const CountMatrix * counts ) const = 0;

};

}

#endif /* REGULARIZOR_H */

