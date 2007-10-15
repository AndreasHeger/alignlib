/*
  alignlib - a library for aligning protein sequences

  $Id: Alignator.h,v 1.3 2004/03/19 18:23:39 aheger Exp $

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

#ifndef ALIGNATOR_H
#define ALIGNATOR_H 1

#include "alignlib.h"

namespace alignlib
{

  class Alignata;
  class Alignandum;
  class Iterator2D;
  class Scorer;
  
  /**  
       @short protocol class for Alignator objects.
       
       Alignator objects align two objects. The default call to use is
       
       Alignator * a;
       Alignata * r;
       
       a->align( row, col, &r );

       The two objects to be aligned to each other are called row and col. This comes from
       the dynamic programming matrix where you have one sequence along the rows and the other 
       along the columns.

       This class is a protocol class and as such defines only the interface.
   
       @author Andreas Heger
       @version $Id: Alignator.h,v 1.3 2004/03/19 18:23:39 aheger Exp $
  */ 
  class Alignator
    {
      /* class member functions-------------------------------------------------------------- */
      
    public:
      /* constructors and desctructors------------------------------------------------------- */
      /** empty constructor */
      Alignator(); 
      
      /** destructor */
      virtual ~Alignator ();
      
      /** copy constructor */
      Alignator( const Alignator & src);
      
      /** method for aligning two arbitrary objects */
      virtual Alignata * align(const Alignandum *, const Alignandum *, Alignata *) = 0;		

      /** accessors */

      /** set range object */
      virtual void setIterator2D( const Iterator2D * iterator = NULL) = 0;    

      /** set scoring object */
      virtual void setScorer( const Scorer * scorer ) = 0;    
      
    };
}

#endif /* ALIGNATOR_H */
