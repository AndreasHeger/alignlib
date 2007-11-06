//--------------------------------------------------------------------------------
// Project LibAlign
//
// Copyright (C) 2000 Andreas Heger All rights reserved
//
// Author: Andreas Heger <heger@ebi.ac.uk>
//
// $Id: ImplTreetorDistanceLinkage.h,v 1.1.1.1 2002/07/08 21:20:17 heger Exp $
//--------------------------------------------------------------------------------

#if HAVE_CONFIG_H
#include <config.h>
#endif

#ifndef IMPL_TREETOR_DISTANCE_LINKAGE_H
#define IMPL_TREETOR_DISTANCE_LINKAGE_H 1

#include <iostream>
#include "alignlib.h"
#include "ImplTreetorDistance.h"

namespace alignlib {
    
/**
   base class for algorithms that generate trees.

   @author Andreas Heger
   @version $Id: ImplTreetorDistanceLinkage.h,v 1.1.1.1 2002/07/08 21:20:17 heger Exp $
   @short asdf
*/ 
class Distor;
class PhyloMatrix;
class ImplTreetorDistanceLinkage : public ImplTreetorDistance {
    
 public:
  /* constructors and desctructors------------------------------------------------------- */

  /** empty constructor */
  ImplTreetorDistanceLinkage ( const Distor * distor, LinkageType method );

  /** empty constructor */
  ImplTreetorDistanceLinkage ( const PhyloMatrix * matrix, LinkageType method );

  /** copy constructor */
  ImplTreetorDistanceLinkage (const ImplTreetorDistanceLinkage & src);
  
  /** destructor */
  virtual ~ImplTreetorDistanceLinkage ();

  /* member access functions--------------------------------------------------------------- */

 protected:
  /** calculate the minimum distance for the working matrix */
  virtual void calculateMinimumDistance() const;

  /** update the distance matrix and other helper variables */
  virtual void updateDistanceMatrix( PhyloMatrixSize cluster_1, PhyloMatrixSize cluster_2 ) const;

  /** join two nodes and return the index of the added node */
  virtual Node joinNodes( PhyloMatrixSize cluster_i, PhyloMatrixSize cluster_2 ) const;

  /** method to use for calculating the distance between new and old clusters */
  LinkageType mMethod;

};

}

#endif /* IMPL_TREETOR_DISTANCE_LINKAGE */

