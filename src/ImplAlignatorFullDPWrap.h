/*
  alignlib - a library for aligning protein sequences

  $Id: ImplAlignatorFullDPWrap.h,v 1.3 2004/03/19 18:23:41 aheger Exp $

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

#ifndef IMPL_ALIGNATOR_FULLDP_WRAP_H
#define IMPL_ALIGNATOR_FULLDP_WRAP_H 1

#include "alignlib.h"
#include "ImplAlignatorFullDP.h"

namespace alignlib {

class SubstitutionMatrix;
class Alignandum;
class Alignata;

/** 
    @short local, full dynamic programming alignment with wrapping around.

    @author Andreas Heger
    @version $Id: ImplAlignatorFullDPWrap.h,v 1.3 2004/03/19 18:23:41 aheger Exp $
*/
class ImplAlignatorFullDPWrap : public ImplAlignatorFullDP {

		      
 public:
    /* Constructors and destructors */

    /** set affine gap penalties and substitution matrix 
     @param subst_matrix	pointer to substitution matrix
     @param row_gop		gap opening penalty in row
     @param row_gep		gap elongation penalty in row
     @param col_gop		gap opening penalty in column, default = row
     @param col_gep		gap elongation penalty in row, default = col
     
    */
    ImplAlignatorFullDPWrap(  const SubstitutionMatrix * subst_matrix ,
			      Score row_gop, Score row_gep, 
			      Score col_gop = 0,Score col_gep = 0 );
			

    /** copy constructor */
    ImplAlignatorFullDPWrap( const ImplAlignatorFullDPWrap & );

    /** destructor */
    virtual ~ImplAlignatorFullDPWrap();

 protected:
    /** perform initialization before alignment */
    virtual void startUp(const Alignandum * row, const Alignandum * col, Alignata * ali);                     
    
    /** perform the alignment */
    virtual void performAlignment(const Alignandum * row, const Alignandum *col, Alignata * result);
};


}

#endif /* _ALIGNATOR_MATRIX_H */

