/*
  alignlib - a library for aligning protein sequences

  $Id: Weightor.h,v 1.3 2004/03/19 18:23:42 aheger Exp $

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

#ifndef WEIGHTOR_H
#define WEIGHTOR_H 1

#include "alignlib.h"

namespace alignlib {

class MultipleAlignment;
class Translator;

/** @short Interface definition for Weightor objects.
    
    Given a multiple alignment of sequences, a Weightor returns a vector of
    weights, for each sequence one weight. The vector has to be deleted
    by the caller.

    These objects take a multiple alignment and return an array of weights.
    Since the same objects will be used for several profiles, no state data 
    between calculations is stored. 
    
    Weighters have to know, what translation was used, because they have to be aware of
    gaps and masked characters. They use the global Translator object for this.
    
    This class is a protocol class and as such defines only the general interface.
    
    @author Andreas Heger
    @version $Id: Weightor.h,v 1.3 2004/03/19 18:23:42 aheger Exp $
*/
class Weightor {
 public:
    // constructors and desctructors

    /** default constructor */
    Weightor();
    
    /** copy constructor */
    Weightor(const Weightor &);

    /** destructor */
    virtual ~Weightor();
    
    /** return a vector of weights for a multiple alignment. The ordering in the result will be the same 
	as in the multiple alignment. Note, that the caller has to delete the weights. */
    virtual SequenceWeight * calculateWeights( const MultipleAlignment & src ) const = 0;

};

}

#endif /* _WEIGHTOR_H */


