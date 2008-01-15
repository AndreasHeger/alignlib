/*
  alignlib - a library for aligning protein sequences

  $Id: ImplRegularizorDirichletFade.h,v 1.2 2004/01/07 14:35:35 aheger Exp $

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

#ifndef IMPL_REGULARIZOR_DIRICHLET_FADE_H
#define IMPL_REGULARIZOR_DIRICHLET_FADE_H 1

#include "alignlib_fwd.h"
#include "ImplRegularizorDirichlet.h"

namespace alignlib {

    typedef double TYPE_WA_COLUMN[PROFILEWIDTH];
    typedef double TYPE_BETA_DIFFERENCES[PROFILEWIDTH][NCOMPONENTS];

    /** When you use this class as a regularizor, the regularizor is 
	switched of, when there are more than FADE_CUTOFF observation 
	in a column.

	@author Andreas Heger
	@version $Id: ImplRegularizorDirichletFade.h,v 1.2 2004/01/07 14:35:35 aheger Exp $
	@short protocol class for sequence weighters
	
    */


class ImplRegularizorDirichletFade : public ImplRegularizorDirichlet {
 public:
    // constructors and desctructors

    /** default constructor */
    ImplRegularizorDirichletFade  ();
    
    /** copy constructor */
    ImplRegularizorDirichletFade  (const ImplRegularizorDirichletFade &);

    /** destructor */
    virtual ~ImplRegularizorDirichletFade ();
    
    /** copy the counts into the frequencies and regularize them by doing so. */
    virtual void fillFrequencies( FrequencyColumn * frequencies, 
				  const CountColumn * counts, 
				  const Position length ) const;

};


}

#endif /* IMPL_REGULARIZOR_DIRICHLET_H */

