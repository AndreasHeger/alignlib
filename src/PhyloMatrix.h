//--------------------------------------------------------------------------------
// Project LibAlign
//
// Copyright (C) 2000 Andreas Heger All rights reserved
//
// Author: Andreas Heger <heger@ebi.ac.uk>
//
// $Id: Matrix.h,v 1.1.1.1 2002/07/08 21:20:17 heger Exp $
//--------------------------------------------------------------------------------

#if HAVE_CONFIG_H
#include <config.h>
#endif

#ifndef PHYLOMATRIX_H
#define PHYLOMATRIX_H 1

#include <iostream>
#include "alignlib.h"

namespace alignlib 
{

/**
   Matrix object, base class for distance and similarity matrices. 
   
   simple and painful
   do not use for heavy computation as the accession operator is
   overloaded. Instead, use matrix template library.
   
   first index is row, second index is column

   @author Andreas Heger
   @version $Id: Matrix.h,v 1.1.1.1 2002/07/08 21:20:17 heger Exp $
   @short contains a matrix
*/ 

class PhyloMatrix 
{
  /* friends---------------------------------------------------------------------------- */
  friend std::ostream & operator<<( std::ostream &, const PhyloMatrix &);
  
  /* class member functions-------------------------------------------------------------- */
 public:
  
  /* constructors and desctructors------------------------------------------------------- */

  /** empty constructor */
  PhyloMatrix ();

  /** copy constructor */
  PhyloMatrix (const PhyloMatrix &);

  /** destructor */
  virtual ~PhyloMatrix ();

  /* member access functions--------------------------------------------------------------- */
  /** return width of PhyloMatrix */
  virtual PhyloMatrixSize getWidth() const = 0;

  /** sets the width of PhyloMatrix, old PhyloMatrix is deleted */
  virtual void setWidth(PhyloMatrixSize width) = 0;

  /** return size (rows * columns) of PhyloMatrix */
  virtual PhyloMatrixSize getSize() const = 0;

  /** return the minimum value of the PhyloMatrix */
  virtual PhyloMatrixValue getMinimum() const = 0;

  /** return the minimum value of the PhyloMatrix + coordinates */
  virtual PhyloMatrixValue getMinimum( Coordinate & coordinates) const = 0;

  /** return the maximum value of the PhyloMatrix */
  virtual PhyloMatrixValue getMaximum() const = 0;

  /** return the maximum value of the PhyloMatrix */
  virtual PhyloMatrixValue getMaximum( Coordinate & coordinates) const = 0;
  
  /** return element */
  virtual PhyloMatrixValue operator()(PhyloMatrixSize row, PhyloMatrixSize col) const = 0;
  virtual PhyloMatrixValue getElement(PhyloMatrixSize row, PhyloMatrixSize col) const = 0;

  /** set element */
  virtual PhyloMatrixValue & operator()(PhyloMatrixSize row, PhyloMatrixSize col) = 0;
  virtual void setElement(PhyloMatrixSize row, PhyloMatrixSize col, PhyloMatrixValue value) = 0;

  /** swap two columns/rows */
  virtual void swap( PhyloMatrixSize col_1, PhyloMatrixSize col_2 ) = 0;

  /** shrink matrix by one */
  virtual void shrink() = 0;
    
  /** read information from stream */
  virtual void read ( std::istream & input ) const = 0;

  virtual void write( std::ostream & output ) const = 0;

};

}

#endif /* _PHYLOMATRIX_H */

