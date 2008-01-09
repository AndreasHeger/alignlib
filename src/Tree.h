//--------------------------------------------------------------------------------
// Project LibAlign
//
// Copyright (C) 2000 Andreas Heger All rights reserved
//
// Author: Andreas Heger <heger@ebi.ac.uk>
//
// $Id: Tree.h,v 1.3 2004/06/02 12:14:35 aheger Exp $
//--------------------------------------------------------------------------------

#if HAVE_CONFIG_H
#include <config.h>
#endif

#ifndef TREE_H
#define TREE_H 1

#include <iosfwd>
#include <string>
#include <vector>
#include "alignlib.h"
#include "alignlib_fwd.h"

/** 
    Base class for trees.

    Trees are graphs, of course. As such, they can be implemented using graph
    libraries (e.g., LEDA, BGL). However, depending on the implementation, 
    different call interfaces are used. Since only one graph library can be
    used at one time, the choice of call interface to the library is used by
    defining WITH_LEDA, WITH_BGL, etc.

    Requirements for a phylogenetic tree:
    1. Store tree structure, i.e. nodes and edges, and allow access methods:
	1. get a Node by Index
	2. get an Edge between two nodes
	3. get Parent of a Node
	4. get left/right Child of a Node
	These methods should be enough to encode custom tree traversal algorithms.

    2. Interact with tree building algorithms. During the calculation of a tree
	the algorithms usually calculate some information that have to be stored
	in the edges/nodes of the tree. These attributes should be generic (using
	templates). 
	Tree building can proceed in three ways:
	1. Bottom up via successive joining of nodes.
	2. Top down via successive addition of nodes to the tree.
	3. Randomly.
	It would be nice to support all three ways of tree construction.

    3. Support tree traversal (needed for Weightors) and updating of node and edge
	properties.

    4. A tree can be either rooted or unrooted. It has to be possible to switch
	between both representations transparently, i.e. so that the traversal 
	algorithms do not need to be aware of it.

    Properties provided by the tree:
    1. num_children contains the number of leaves below a node. The num_children
	of a leaf is 1, of the root it is equal to getNumLeaves(). This property
	is set automatically when joinNodes is called.

    2. height measures the height of a node above the leaves. This property has to
	be set by the user, since there are different ways to calculate the height
	of a node.

    This class is a protocol class and as such defines only the general interface.

    @author Andreas Heger
    @version $Id: Tree.h,v 1.3 2004/06/02 12:14:35 aheger Exp $
    @short protocol class for phylogenetic trees.
*/
 
namespace alignlib 
{

class Tree 
{
  
  /* friends---------------------------------------------------------------------------- */
  friend std::ostream & operator<<( std::ostream &, const Tree &);

  /* class member functions-------------------------------------------------------------- */
 public:
  
  /* constructors and desctructors------------------------------------------------------- */
  Tree();

  /** copy constructor */
  Tree (const Tree & src);
  
  /** destructor */
  virtual ~Tree ();

  //------------------------------------------------------------------------------------------------------------
  /** return a new object of the same type */
  virtual HTree getNew() const = 0;
  
  /** return an identical copy */
  virtual HTree getClone() const = 0;
  
  /* member access functions--------------------------------------------------------------- */

  /** returns the number of leaves */
  virtual Node getNumLeaves() const = 0;

  /** sets the number of leaves. This erases the Tree and allocates memory 
      for a new one with num_leaves leaves */  
  virtual void setNumLeaves( unsigned int num_leaves ) = 0;

  /** set the height of a vertex */
  virtual void setHeight( Node node, TreeHeight height ) = 0;

  /** get the height of a vertex */
  virtual TreeHeight getHeight( Node node ) const = 0;

  /** set the weight of an edge */
  virtual void setWeight( Node child, Node parent, TreeWeight weight ) = 0;

  /** get the weight of an edge */
  virtual TreeWeight getWeight( Node child, Node parent ) const = 0;

  /** returns the root */
  virtual Node getNoNode() const = 0;

  /** returns the root */
  virtual Node getRoot() const = 0;

  /** returns the parent of a node */
  virtual Node getParent( Node node ) const = 0;

  /** returns the left child of a node */
  virtual Node getLeftChild( Node node ) const = 0;

  /** returns the right of a node */
  virtual Node getRightChild( Node node ) const = 0;

  /** returns the number of children of a node */
  virtual Node getNumChildren( Node node ) const = 0;

  /** returns a vector of leaves nodes */
  virtual HNodeVector getNodesLeaves() const = 0;

  /** returns a vector of nodes sorted according to breadth-first-traversal, first encounter */
  virtual HNodeVector getNodesBreadthFirstFinish() const = 0;

  /** returns a vector of nodes sorted according to breadth-first-traversal, last encounter */
  virtual HNodeVector getNodesBreadthFirstVisit() const = 0;

  /** returns a vector of nodes sorted according to depth-first-traversal, first encounter */
  virtual HNodeVector getNodesDepthFirstVisit() const = 0;

  /** returns a vector of nodes sorted according to depth-first-traversal, last encounter */
  virtual HNodeVector getNodesDepthFirstFinish() const = 0;

  /* ---------------------------------------------------------------------------------------- */
  /** removes the root */
  virtual void removeRoot() = 0 ;

  /** sets the root */
  virtual Node setRoot( const Node node_1, 
			     const Node node_2,
			     TreeWeight weight ) = 0;

  /** find root */
  virtual Node findRoot( const Node node ) const = 0;

  /** Add a node to a tree for bottom-up construction, i.e. a construction which starts by
      joining leaves. This function returns a reference to the internal node that was just 
      created from two already existing nodes.

      @param node_1	node which is joined
      @param node_2	node which is joined
      @param weight_1	weight of edge 1 to new internal node
      @param weight_2	weight of edge 2 to new internal node
  */
  virtual Node joinNodes( const Node node_1, 
			       const Node node_2,
			       TreeWeight weight_1, 
			       TreeWeight weight_2,
			       const bool map_parents = false ) = 0;

  virtual void write( std::ostream & output ) const = 0;

};

}

#endif /* TREE_H */

