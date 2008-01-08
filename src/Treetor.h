//--------------------------------------------------------------------------------
// Project LibAlign
//
// Copyright (C) 2000 Andreas Heger All rights reserved
//
// Author: Andreas Heger <heger@ebi.ac.uk>
//
// $Id: Treetor.h,v 1.1.1.1 2002/07/08 21:20:17 heger Exp $
//--------------------------------------------------------------------------------

#if HAVE_CONFIG_H
#include <config.h>
#endif

#ifndef TREETOR_H
#define TREETOR_H 1

namespace alignlib 
{
    class MultipleAlignment;

/**
   base class for algorithms that generate trees. So far this class is empty

   This class is a protocoll class and as such only defines an empty
   interface.

   @author Andreas Heger
   @version $Id: Treetor.h,v 1.1.1.1 2002/07/08 21:20:17 heger Exp $
   @short Algorithm class that generates trees.
*/ 

class Tree;

class Treetor {

  /* class member functions-------------------------------------------------------------- */
 public:

  /* constructors and desctructors------------------------------------------------------- */
  /** empty constructor */
  Treetor ();

  /** copy constructor */
  Treetor (const Treetor & src);
  
  /** destructor */
  virtual ~Treetor ();

  /** create a tree from a multiple alignment */
  virtual Tree * calculateTree( Tree * dest, const alignlib::HMultipleAlignment src = NULL) const = 0; 

};

}

#endif /* TREETOR_H */

