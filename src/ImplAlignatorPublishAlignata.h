/*
  alignlib - a library for aligning protein sequences

  $Id: ImplAlignatorPublishAlignata.h,v 1.2 2004/01/07 14:35:34 aheger Exp $

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

#ifndef IMPL_ALIGNATOR_PUBLISHALIGNATA_H
#define IMPL_ALIGNATOR_PUBLISHALIGNATA_H 1

#include "alignlib.h"
#include "ImplAlignator.h"

namespace alignlib {

 /** This alignator object is initialized with an alignment. 
     The alignment is returned after calling the method align.
     Use this for storing premade dot-plots.
     
     @author Andreas Heger
     @version $Id: ImplAlignatorPublishAlignata.h,v 1.2 2004/01/07 14:35:34 aheger Exp $
     @short Return a non-const pre-built alignment
      
  */

class Alignandum;

class ImplAlignatorPublishAlignata : public ImplAlignator 
{
 public:
    // constructors and desctructors

    /** default constructor */
    ImplAlignatorPublishAlignata  ( Alignata * ali);
    
    /** copy constructor */
    ImplAlignatorPublishAlignata  (const ImplAlignatorPublishAlignata &);

    /** destructor */
    virtual ~ImplAlignatorPublishAlignata();

    /** method for aligning two arbitrary objects */
    virtual Alignata * align(const Alignandum *, const Alignandum *, Alignata *);

    /** return a new alignator object of the same type.
     */
    virtual ImplAlignatorPublishAlignata * getClone() const;
    
 private:
    Alignata * mAlignata;

};

}

#endif /* IMPL_DOTTOR_H */

