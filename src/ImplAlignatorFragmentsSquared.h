/*
  alignlib - a library for aligning protein sequences

  $Id: ImplAlignatorFragmentsSquared.h,v 1.3 2004/03/19 18:23:41 aheger Exp $

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

#ifndef IMPL_ALIGNATOR_FRAGMENTS_SQUARED_H
#define IMPL_ALIGNATOR_FRAGMENTS_SQUARED_H 1

#include "alignlib.h"
#include "ImplAlignatorFragments.h"

namespace alignlib {
    class Alignandum;
    class Alignata;
    class SubstitutionMatrix;
    class Fragmentor;
    
/** 
    This objects aligns two @ref Alignandum objects using a dot-matrix.
    
    @author Andreas Heger
    @version $Id: ImplAlignatorFragmentsSquared.h,v 1.3 2004/03/19 18:23:41 aheger Exp $
    @short represents an alignment algorithm
*/
class ImplAlignatorFragmentsSquared : public ImplAlignatorFragments {
 public:

    /* constructors and destructors */

    /** set affine gap penalties
     @param row_gop		gap opening penalty in row
     @param row_gep		gap elongation penalty in row
     @param col_gop		gap opening penalty in column, default = row
     @param col_gep		gap elongation penalty in row, default = col
     
    */
    ImplAlignatorFragmentsSquared( Score row_gop, 
				   Score row_gep, 
				   Score col_gop = 0,
				   Score col_gep = 0,
				   Fragmentor * fragmentor = NULL );

    /** copy constructor */
    ImplAlignatorFragmentsSquared( const ImplAlignatorFragmentsSquared & );

    /** destructor */
    virtual ~ImplAlignatorFragmentsSquared();
    
    /** return a new alignator object of the same type.
     */
    virtual ImplAlignatorFragmentsSquared * getClone() const;    

 protected:

    /** perform the alignment */
    virtual void performAlignment(const Alignandum * row, const Alignandum *col, Alignata * result);

    /** get GAP cost for a gap */
    virtual Score getGapCost( Dot x1, Dot x2 ) const;
};


}

#endif /* _ALIGNATOR_DOTS_H */

