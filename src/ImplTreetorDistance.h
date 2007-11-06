//--------------------------------------------------------------------------------
// Project LibAlign
//
// Copyright (C) 2000 Andreas Heger All rights reserved
//
// Author: Andreas Heger <heger@ebi.ac.uk>
//
// $Id: ImplTreetorDistance.h,v 1.1.1.1 2002/07/08 21:20:17 heger Exp $
//--------------------------------------------------------------------------------

#if HAVE_CONFIG_H
#include <config.h>
#endif

#ifndef IMPL_TREETOR_DISTANCE_H
#define IMPL_TREETOR_DISTANCE_H 1

#include <iosfwd>

#include "alignlib.h"
#include "ImplTreetor.h"

namespace alignlib {
  class MultipleAlignment;
}

namespace alignlib {
    
/**
   Base class for algorithms that generate trees based on distance matrices.
   
   This class needs a distance matrix to work on. The matrix is copied, when
   it is supplied in the factory function.

   @author Andreas Heger
   @version $Id: ImplTreetorDistance.h,v 1.1.1.1 2002/07/08 21:20:17 heger Exp $
   @short contains a matrix
*/ 

class Tree;
class Distor;
class PhyloMatrix;

class ImplTreetorDistance : public ImplTreetor 
{

  /* class member functions-------------------------------------------------------------- */
 public:
  /* constructors and desctructors------------------------------------------------------- */

  /** constructor */
  ImplTreetorDistance( const Distor * distor );

  /** copy constructor */
  ImplTreetorDistance (const ImplTreetorDistance & src);
  
  /** destructor */
  virtual ~ImplTreetorDistance ();

  /* member access functions--------------------------------------------------------------- */
  /** create a tree with distance in matrix. The algorithm is generic, the different clustering
   strategies are implemented using the functions CalculateMinimunDistance, Initialize, Cleanup, ...
   @param tree	pointer to tree object, where result is stored. If it is not supplied, the method
   creates one and returns it. This can be useful to pass different tree implementations. However the
   current implementation is optimal for building trees from distance matrices.
  */

  /** create a tree from a multiple alignment */
  virtual Tree * calculateTree( Tree * dest, const alignlib::MultipleAlignment * src) const; 

 protected:
  /** initialize helper variables that you might need */
  virtual void startUp(Tree * dest, const alignlib::MultipleAlignment * src) const;

  /** clean up helper variables */
  virtual void cleanUp() const;

  /** swap indices in helper arrays */
  virtual void swapHelpers( PhyloMatrixSize cluster_1, PhyloMatrixSize cluster_2) const;

  /** calculate the minimum distance for the working matrix */
  virtual void calculateMinimumDistance() const = 0;

  /** update the distance matrix and other helper variables */
  virtual void updateDistanceMatrix( PhyloMatrixSize cluster_1, PhyloMatrixSize cluster_2 ) const = 0;

  /** join two nodes and return the index of the added node */
  virtual Node joinNodes( PhyloMatrixSize cluster_i, PhyloMatrixSize cluster_2 ) const = 0;

    // member data
 public:

 protected:

  /** Distor object used for calculating distance matrix */
  const Distor * mDistor;

  /** coordinates of minimum element in distance matrix
   */
  mutable Coordinate mMinimumCoord;

  /** minimum distance in distance matrix
   */
  mutable PhyloMatrixValue mMinimumValue;

  /** working copy of distance matrix */
  mutable PhyloMatrix * mWorkMatrix;
  
  /** indices, mapping rows in distance matrix to nodes in tree*/
  mutable Node * mIndices;

  /** pointer to result, used for internal access between methods */
  mutable Tree * mTree;
  
};

}

#endif /* IMPL_TREETOR_DISTANCE_H */

