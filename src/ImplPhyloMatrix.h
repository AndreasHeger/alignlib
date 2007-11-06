//--------------------------------------------------------------------------------
// Project LibAlign
//
// Copyright (C) 2000 Andreas Heger All rights reserved
//
// Author: Andreas Heger <heger@ebi.ac.uk>
//
// $Id: ImplMatrix.h,v 1.1.1.1 2002/07/08 21:20:17 heger Exp $
//--------------------------------------------------------------------------------

#if HAVE_CONFIG_H
#include <config.h>
#endif

#ifndef IMPL_PHYLOMATRIX_H
#define IMPL_PHYLOMATRIX_H 1

#include <iosfwd>

#include "PhyloMatrix.h"

namespace alignlib {

/**

   -> memory is allocated when the matrix is created and not given free, until it
   is deleted (i.e. no partial memory freeing during shrinks)

   @author Andreas Heger
   @version $Id: ImplMatrix.h,v 1.1.1.1 2002/07/08 21:20:17 heger Exp $
   @short contains a matrix
*/ 

  // i.e., for invalid requests return first element. This is bad, but
  // did not want to do all the checking
#define NO_INDEX 0		

class ImplPhyloMatrix : public PhyloMatrix {
 public:
  
  /* constructors and desctructors------------------------------------------------------- */

  /** empty constructor */
  ImplPhyloMatrix ();

  /** create PhyloMatrix of width and set all values to default_value */
  ImplPhyloMatrix( PhyloMatrixSize width, PhyloMatrixValue default_value);

  /** load array */
  ImplPhyloMatrix( PhyloMatrixSize width, PhyloMatrixValue * source);

  /** copy constructor */
  ImplPhyloMatrix (const ImplPhyloMatrix &);

  /** destructor */
  virtual ~ImplPhyloMatrix ();

  /* member access functions--------------------------------------------------------------- */
  /** return width of PhyloMatrix */
  virtual PhyloMatrixSize getWidth() const;

  /** sets the width of PhyloMatrix, old PhyloMatrix is deleted */
  virtual void setWidth(PhyloMatrixSize width);

  /** return size (rows * columns) of PhyloMatrix */
  virtual PhyloMatrixSize getSize() const;

  /** return the minimum value of the PhyloMatrix */
  virtual PhyloMatrixValue getMinimum() const;

  /** return the minimum value of the PhyloMatrix + coordinates */
  virtual PhyloMatrixValue getMinimum( Coordinate & coordinates) const;

  /** return the maximum value of the PhyloMatrix */
  virtual PhyloMatrixValue getMaximum() const;

  /** return the maximum value of the PhyloMatrix */
  virtual PhyloMatrixValue getMaximum( Coordinate & coordinates) const;
  
  /** return element */
  virtual PhyloMatrixValue operator()(PhyloMatrixSize row, PhyloMatrixSize col) const;
  virtual PhyloMatrixValue getElement(PhyloMatrixSize row, PhyloMatrixSize col) const;

  /** set element */
  virtual PhyloMatrixValue & operator()(PhyloMatrixSize row, PhyloMatrixSize col);
  virtual void setElement(PhyloMatrixSize row, PhyloMatrixSize col, PhyloMatrixValue value);

  /** delete last row/column from PhyloMatrix */
  virtual void Shrink();

  /** swap two columns/rows */
  virtual void Swap( PhyloMatrixSize col_1, PhyloMatrixSize col_2 );

  /** read information from stream */
  virtual void Read ( std::istream & input ) const;

  virtual void Write( std::ostream & output ) const;


 protected:
  
  /** get row for given index */
  virtual PhyloMatrixSize getRow( PhyloMatrixSize index ) const;
 
  /** get column for given index */
  virtual PhyloMatrixSize getColumn( PhyloMatrixSize index ) const;
 
  /** return mapped index for row and column */
  virtual PhyloMatrixSize getIndex(PhyloMatrixSize row, PhyloMatrixSize col) const;   

  /** allocate memory */
  virtual void allocateMemory();

  /** allocate memory */
  virtual void freeMemory();

  /** calculate the memory size of the PhyloMatrix given its width*/
  virtual void calculateSize();

  /** the width of the PhyloMatrix */
  PhyloMatrixSize mWidth;
 
  /** the size of the PhyloMatrix in terms of elements of type TYPE_PhyloMatrix*/
  PhyloMatrixSize mSize;

  /** a pointer to the matrix */
  PhyloMatrixValue * mMatrix;

};


}

#endif /* _ALIGNATUMSEQUENCE_H */

