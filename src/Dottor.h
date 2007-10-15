/*
  alignlib - a library for aligning protein sequences

  $Id: Dottor.h,v 1.3 2004/03/19 18:23:40 aheger Exp $

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

#ifndef DOTTOR_H
#define DOTTOR_H 1

#include "alignlib.h"

namespace alignlib 
{
  
class Alignandum;
class SubstitutionMatrix;
struct ResiduePAIR;

/** Protocoll class for objects that calculate aligned residue pairs (dots) for two
    sequences
    
    This class is a protocol class and as such defines only the general interface.

    This class is a helper class for Alignator objects, that need a list of
    dots for calculating an alignment.
    
    @author Andreas Heger
    @version $Id: Dottor.h,v 1.3 2004/03/19 18:23:40 aheger Exp $
    @short protocol class for functions that calculate dots (aligned residue pairs)
    
*/
class Dottor {
 public:
    // constructors and desctructors

    /** default constructor */
    Dottor  ();
    
    /** copy constructor */
    Dottor  (const Dottor &);

    /** destructor */
    virtual ~Dottor ();

    /** release memory hold by Dottor */
    virtual void releasePairs() const = 0;
    
    /** calculate the dots */
    virtual void calculatePairs( const Alignandum * row, const Alignandum * col) const = 0;
    
    /** get array of residue pairs */
    virtual ResiduePAIR * getPairs() const = 0;

    /** get location of array for row-indices */
    virtual Position * getRowIndices() const = 0;

};


/** create a DottorIdentity object.

    This object will calculate dots based on residue identity.
 */
Dottor * makeDottorIdentity();

 
/** create a DottorSimilarity object.

    This object will calculate dots based on pairwise similarity.
    
    @param matrix    Substitution matrix for determining similarity.
*/
Dottor * makeDottorSimilarity( const SubstitutionMatrix * matrix);
 

/** gets the default Dottor object */ 
const Dottor * getDefaultDottor();

/** sets the default Dottor object */
const Dottor * setDefaultDottor( const Dottor * dottor);


}

#endif /* DOTTOR_H */

