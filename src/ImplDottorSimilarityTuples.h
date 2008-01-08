/*
  alignlib - a library for aligning protein sequences

  $Id: ImplDottorSimilarityTuples.h,v 1.2 2004/01/07 14:35:35 aheger Exp $

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

#ifndef IMPL_DOTTOR_SIMILARITY_TUPLES_H
#define IMPL_DOTTOR_SIMILARITY_TUPLES_H 1

#include "alignlib.h"
#include "ImplDottorSimilarity.h"

namespace alignlib {

  /** Protocoll class for objects, that regularize profile columns.
      
      Create a dot for similar residue pairs. Similarity is given by
      a substitution matrix.
      
      @author Andreas Heger
      @version $Id: ImplDottorSimilarityTuples.h,v 1.2 2004/01/07 14:35:35 aheger Exp $
      @short Create a dot-matrix using identical residue pairs.
      
  */

class Alignandum;
class SubstitutionMatrix;

class ImplDottorSimilarityTuples : public ImplDottor {
 public:
    // constructors and desctructors

    /** default constructor */
    ImplDottorSimilarityTuples  (const SubstitutionMatrix * matrix, int ktuple);
    
    /** copy constructor */
    ImplDottorSimilarityTuples  (const ImplDottorSimilarityTuples &);

    /** destructor */
    virtual ~ImplDottorSimilarityTuples ();

 protected:

    /** calculate the dots. Overload this method to calculate different types of dots */
    virtual void calculateNewPairs( const HAlignandum row, const HAlignandum col) const;
    
};

}

#endif /* IMPL_DOTTOR_H */

