/*
  alignlib - a library for aligning protein sequences

  $Id: ImplFragmentorDiagonals.h,v 1.3 2004/03/19 18:23:41 aheger Exp $

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

#ifndef IMPL_FRAGMENTOR_DIAGONALS_H
#define IMPL_FRAGMENTOR_DIAGONALS_H 1

#include "alignlib.h"
#include "ImplFragmentor.h"
#include "Alignment.h"
#include "SubstitutionMatrix.h"

namespace alignlib 
{

	class Alignandum;
	class Alignment;

/**
   @short align fragments using a gap penalty based on diagonal changes.
   
   @author Andreas Heger
   @version $Id: ImplFragmentorDiagonals.h,v 1.3 2004/03/19 18:23:41 aheger Exp $

*/ 
class ImplFragmentorDiagonals : public ImplFragmentor {
  /* class member functions-------------------------------------------------------------- */
 public:
    /* constructors and desctructors------------------------------------------------------- */
    /** constructor */
    ImplFragmentorDiagonals( Score row_gop, 
			     Score row_gep, 
			     Score col_gop = 0 , 
			     Score col_gep = 0,
			     Alignator * dottor = NULL);

    /** destructor */
    virtual ~ImplFragmentorDiagonals ();

    /** copy constructor */
    ImplFragmentorDiagonals( const ImplFragmentorDiagonals & src);

 protected:
    /** gap opening penalty for row-object */
    Score mRowGop;
    /** gap elongation penalty for col-object */
    Score mRowGep;
    /** gap opening penalty for row-object */
    Score mColGop;
    /** gap elongation penalty for col-object */
    Score mColGep;  

    /** the alignator used to create dot-plots */
    Alignator *mDottor;

    virtual Score getGapCost( const ResiduePAIR & p1, const ResiduePAIR & p2 ) const ;

    /** perform the actual alignment */
    virtual void performFragmentation( const Alignandum * row, 
				       const Alignandum * col, 
				       const Alignment * sample);
    
    /** perform cleanup after alignment */
    virtual void cleanUp(const Alignandum * row, const Alignandum * col, Alignment * ali);                     

};

}

#endif /* FRAGMENTOR_H */
