//--------------------------------------------------------------------------------
// Project LibAlign
//
// Copyright (C) 2000 Andreas Heger All rights reserved
//
// Author: Andreas Heger <heger@ebi.ac.uk>
//
// $Id: ImplPhyloMatrixSymmetric.h,v 1.1.1.1 2002/07/08 21:20:17 heger Exp $
//--------------------------------------------------------------------------------

#if HAVE_CONFIG_H
#include <config.h>
#endif

#ifndef IMPL_PHYLOMATRIX_SYMMETRIC_H
#define IMPL_PHYLOMATRIX_SYMMETRIC_H 1

#include <iosfwd>

#include "ImplPhyloMatrix.h"

namespace alignlib {

/**

   @author Andreas Heger
   @version $Id: ImplPhyloMatrixSymmetric.h,v 1.1.1.1 2002/07/08 21:20:17 heger Exp $
   @short contains a PhyloMatrix
*/ 

class ImplPhyloMatrixSymmetric : public ImplPhyloMatrix {
 public:
  
  /* constructors and desctructors------------------------------------------------------- */

  /** empty constructor */
  ImplPhyloMatrixSymmetric ();

  /** create PhyloMatrix of width and set all values to default_value */
  ImplPhyloMatrixSymmetric( PhyloMatrixSize width, PhyloMatrixValue default_value);

  /** load array */
  ImplPhyloMatrixSymmetric( PhyloMatrixSize width, PhyloMatrixValue * source);

  /** copy constructor */
  ImplPhyloMatrixSymmetric (const ImplPhyloMatrixSymmetric &);

  /** destructor */
  virtual ~ImplPhyloMatrixSymmetric ();

  /* member access functions--------------------------------------------------------------- */

  /** delete last row/column from PhyloMatrix */
  virtual void Shrink();

  /** swap two columns/rows */
  virtual void Swap( PhyloMatrixSize col_1, PhyloMatrixSize col_2 );

 protected:
  
  /** get row for given index */
  virtual PhyloMatrixSize getRow( PhyloMatrixSize index ) const;
 
  /** get column for given index */
  virtual PhyloMatrixSize getColumn( PhyloMatrixSize index ) const;
 
  /** return mapped index for row and column */
  virtual PhyloMatrixSize getIndex(PhyloMatrixSize row, PhyloMatrixSize col) const;   

  /** calculate the size of the PhyloMatrix */
  virtual void calculateSize();

};


}

#endif /* _ALIGNATUMSEQUENCE_H */

