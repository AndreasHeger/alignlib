/*
  alignlib - a library for aligning protein sequences

  $Id: ImplWeightorHenikoff.h,v 1.3 2004/03/19 18:23:41 aheger Exp $

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

#ifndef IMPL_WEIGHTOR_HENIKOFF_H
#define IMPL_WEIGHTOR_HENIKOFF_H 1

#include "alignlib.h"
#include "ImplWeightor.h"

namespace alignlib {

class MultipleAlignment;
class Translator;

/** @short Implements sequence weighting as proposed by the Henikoffs.
    
    This weighter uses the weighting scheme by Henikoff. Weights are normalized, so 
    that they sum to 1.
     
    @author Andreas Heger
    @version $Id: ImplWeightorHenikoff.h,v 1.3 2004/03/19 18:23:41 aheger Exp $
*/
class ImplWeightorHenikoff : public ImplWeightor {
 public:
    // constructors and desctructors

    /** default constructor */
    ImplWeightorHenikoff(const Translator * translator);
    
    /** copy constructor */
    ImplWeightorHenikoff(const ImplWeightorHenikoff &);

    /** destructor */
    virtual ~ImplWeightorHenikoff();
    
    /** return a vector of weights for a multiple alignment. The ordering in the result will be the same 
	as in the multiple alignment. Note, that the caller has to delete the weights. */
    virtual SequenceWeight * calculateWeights( const MultipleAlignment & src ) const;

};


}

#endif /* IMPL_WEIGHTOR_HENIKOFF_H */

