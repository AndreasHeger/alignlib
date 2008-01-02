/*
  alignlib - a library for aligning protein sequences

  $Id: ImplFragmentorRepetitive.h,v 1.3 2004/03/19 18:23:41 aheger Exp $

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

#ifndef IMPL_FRAGMENTOR_REPETITIVE_H
#define IMPL_FRAGMENTOR_REPETITIVE_H 1

#include "alignlib.h"
#include "ImplFragmentor.h"
#include "Alignment.h"
#include "SubstitutionMatrix.h"

namespace alignlib 
{

	class Alignandum;
	class Alignment;
	class Alignator;

/**
   @short build fragments by repeatedly applying an Alignator object.

   Previously aligned regions are masked.
   
   @author Andreas Heger
   @version $Id: ImplFragmentorRepetitive.h,v 1.3 2004/03/19 18:23:41 aheger Exp $
*/ 

class ImplFragmentorRepetitive : public ImplFragmentor {
  /* class member functions-------------------------------------------------------------- */
 public:
    /* constructors and desctructors------------------------------------------------------- */
    /** constructor */
    ImplFragmentorRepetitive( Alignator * alignator, 
			      Score min_score );

    /** destructor */
    virtual ~ImplFragmentorRepetitive ();

    /** copy constructor */
    ImplFragmentorRepetitive( const ImplFragmentorRepetitive & src);

 protected:
    /** the alignator used to create dot-plots */
    Alignator * mAlignator;

    /** minimum score for alignment */
    Score mMinScore;

    /** perform the actual alignment */
    virtual void performFragmentation( const Alignandum * row, 
				       const Alignandum * col, 
				       const Alignment * sample);
};

}

#endif /* FRAGMENTOR_H */
