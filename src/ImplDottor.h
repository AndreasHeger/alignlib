/*
  alignlib - a library for aligning protein sequences

  $Id: ImplDottor.h,v 1.3 2004/03/19 18:23:41 aheger Exp $

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

#ifndef IMPL_DOTTOR_H
#define IMPL_DOTTOR_H 1

#include "alignlib_fwd.h"
#include "Dottor.h"

namespace alignlib {

class Alignandum;

class ImplAlignmentMatrix; 

  /** Protocoll class for objects, that regularize profile columns.
      
      This class is a protocol class and as such defines only the general interface.
      
      @author Andreas Heger
      @version $Id: ImplDottor.h,v 1.3 2004/03/19 18:23:41 aheger Exp $
      @short protocol class for sequence weighters
      
  */
class ImplDottor : public Dottor {
 public:
    // constructors and desctructors

    /** default constructor */
    ImplDottor  ();
    
    /** copy constructor */
    ImplDottor  (const ImplDottor &);

    /** destructor */
    virtual ~ImplDottor ();

    /** release memory hold by ImplDottor */
    virtual void releasePairs() const;
    
    /** calculate the dots */
    virtual void calculatePairs( const HAlignandum row, const HAlignandum col) const;
    
    /** get array of residue pairs */
    virtual ResiduePAIR * getPairs() const;

    /** get location of array for row-indices */
    virtual Position * getRowIndices() const;

 protected:
    mutable ImplAlignmentMatrix * mMatrix;

    /** calculate the dots. Overload this method to calculate different types of dots */
    virtual void calculateNewPairs( const HAlignandum row, const HAlignandum col) const = 0;

};

}

#endif /* IMPL_DOTTOR_H */

