//--------------------------------------------------------------------------------
// Project alignlib
//
// Copyright (C) 2000 Andreas Heger All rights reserved
//
// Author: Andreas Heger <heger@ebi.ac.uk>
//
// $Id: ImplTree.cpp,v 1.3 2004/06/02 12:14:34 aheger Exp $
//--------------------------------------------------------------------------------

#include <iostream>
#include <iomanip>
#include <string>
#include <utility>
#include <stack>
#include <cassert>
#include "ImplTree.h"
#include "AlignlibDebug.h"

using namespace std;

namespace alignlib 
{

//-----------------------------< factory functions >----------------------------------------------------
Tree * makeTree( Node num_leaves) { return new ImplTree( num_leaves); }

std::ostream & operator<<( std::ostream & output, const NODE_INFO & src) {
  output << src.mParent << " " << src.mLeftChild << " " << src.mRightChild << " " 
	 << src.mNumChildren << " " << src.mWeight << " " << src.mHeight;

  return output;
}

//---------------------------------------------------------< constructors and destructors >--------------------------------------
ImplTree::ImplTree () : mNumLeaves(0), mTree(NULL) {
}

ImplTree::ImplTree( size_t num_leaves) : 
  mNumLeaves( num_leaves ), mCurrentNode(0), mTree(NULL) {
  
  if (num_leaves > 0) 
    setNumLeaves( num_leaves );
}
		       
ImplTree::~ImplTree () 
{
	debug_func_cerr(5);
  delete [] mTree;
}

ImplTree::ImplTree (const ImplTree & src ) : 
  mNumLeaves ( src.mNumLeaves ), mCurrentNode( src.mCurrentNode ) 
  {
	debug_func_cerr(5);

	if (mTree != NULL) delete [] mTree;
	
	if (src.mTree != NULL)
	{
		size_t n = mNumLeaves * 2 - 1;
		mTree = new NODE_INFO[ n ];
		memcpy( mTree, src.mTree, sizeof(NODE_INFO) * n );
	}
		
}

//-------------------------------------------------------< accessors >----------------------------------------
Node ImplTree::getNumLeaves() const {
  return mNumLeaves;
}    

//-------------------------------------------------------< accessors >----------------------------------------
void ImplTree::setNumLeaves(unsigned int num_leaves) {

	debug_func_cerr(5);

  assert(num_leaves != 0);
  delete [] mTree;
  mNumLeaves = num_leaves;
  mTree = new NODE_INFO[num_leaves * 2 - 1];
  mCurrentNode = mNumLeaves;
  recordLeaves();
}    

Node ImplTree::getNoNode()  const {
  return NO_NODE;
}
  
Node ImplTree::getRoot()  const {
  return mCurrentNode - 1;
}

Node ImplTree::getLeftChild( Node node) const {
  return mTree[node].mLeftChild;
}

Node ImplTree::getRightChild( Node node) const {
  return mTree[node].mRightChild;
}

Node ImplTree::getParent( Node node) const {
  return mTree[node].mParent;
}

Node ImplTree::getNumChildren( Node node) const {
  return mTree[node].mNumChildren;
}

TreeHeight ImplTree::getHeight( Node node) const {
  return mTree[node].mHeight;
}

void ImplTree::setHeight( Node node, TreeHeight height) {
  mTree[node].mHeight = height;
}

TreeWeight ImplTree::getWeight( Node child, Node parent) const 
{
    return mTree[child].mWeight;
}

void ImplTree::setWeight( Node child, Node parent, TreeWeight weight) 
{
     mTree[child].mWeight = weight;
}

//-------------------------------------------------------< others >-------------------------------------------
Node ImplTree::findRoot( const Node node) const 
{
  
  Node n = node;
  while( mTree[n].mParent != NO_NODE)
    n = mTree[n].mParent;
  return n;
}

Node ImplTree::joinNodes( const Node node_1, 
			       const Node node_2,
			       const TreeWeight edge_weight_1, 
			       const TreeWeight edge_weight_2,
			       const bool map_parents) 
{

	debug_func_cerr( 5 );
	
    Node node1, node2;

    if (map_parents) {
      node1 = findRoot( node_1 );
      node2 = findRoot( node_2 );
    } else {
      node1 = node_1;
      node2 = node_2;
    }
    
    mTree[node1].mParent = mCurrentNode;
    mTree[node2].mParent = mCurrentNode;
    
    mTree[node1].mWeight = edge_weight_1;
    mTree[node2].mWeight = edge_weight_2;

    mTree[mCurrentNode].mLeftChild = node1;
    mTree[mCurrentNode].mRightChild = node2;
    mTree[mCurrentNode].mNumChildren = mTree[node1].mNumChildren + 
      mTree[node2].mNumChildren;
    
    return mCurrentNode++;
}

//------------------------------------------------------------------------------------------------------------
void ImplTree::removeRoot() 
{
	debug_func_cerr( 5 );
}
  
//------------------------------------------------------------------------------------------------------------
Node ImplTree::setRoot( Node node_from, Node node_to, TreeWeight weight ) 
{
	debug_func_cerr( 5 );
	return getRoot();
}

//----------------------------------------------------------------------------------------
void ImplTree::recordLeaves() 
{
	debug_func_cerr( 5 );

  // leaves are the first nodes in the graph starting at index 0
  for (Node node = 0; node < getNumLeaves(); node ++) 
  {
    
    mTree[node].mHeight = 0;
    mTree[node].mWeight = 0;
    mTree[node].mNumChildren = 1;
    mTree[node].mParent = NO_NODE;
    mTree[node].mLeftChild = NO_NODE;
    mTree[node].mRightChild = NO_NODE;
  }
}


//---------------------------------> methods for creation of node-lists <--------------------------------------------
//--------------------------------------------------------------------------------------------------------------------
NodeVector * ImplTree::getNodesLeaves() const 
{
	debug_func_cerr( 5 );

    std::vector<Node> * nodes = new std::vector<Node>(getNumLeaves());

    // by definition, all nodes in the graph with number 0 to getNumLeaves - 1 are leaves
    // there is probably a STL version of this
    for (Node i = 0; i < getNumLeaves(); i++)
	(*nodes)[i] = i;

    return nodes;
}
    
//--------------------------------------------------------------------------------------------------------------------
NodeVector * ImplTree::getNodesDepthFirstVisit() const 
{
	debug_func_cerr( 5 );

  std::vector<Node> * nodes = new std::vector<Node>;

  std::stack< Node> s;
  
  s.push( getRoot() );
  
  while (!s.empty()) {
    
    Node node = s.top();
    s.pop();
    
    nodes->push_back( node );
    
    Node c;
    if ((c = getLeftChild(node)) != NO_NODE)
      s.push( c );
    if ((c = getRightChild(node)) != NO_NODE)
      s.push( c );
  }
  
  return nodes;
}

//--------------------------------------------------------------------------------------------------------------------
void ImplTree::traversePostOrder( Node node, NodeVector * nodes)  const 
{
	debug_func_cerr( 5 );
  
  if (node == NO_NODE)
    return;

  traversePostOrder(getLeftChild(node), nodes);
  traversePostOrder(getRightChild(node), nodes);

  nodes->push_back( node );
}


NodeVector * ImplTree::getNodesDepthFirstFinish() const 
{
	debug_func_cerr( 5 );

  std::vector<Node> * nodes = new std::vector<Node>;

  traversePostOrder( getRoot(), nodes);
  return nodes;

};

//--------------------------------------------------------------------------------------------------------------------
NodeVector * ImplTree::getNodesBreadthFirstVisit() const 
{
	debug_func_cerr( 5 );
	
  assert( 0 );
  return NULL;
}

//--------------------------------------------------------------------------------------------------------------------
NodeVector * ImplTree::getNodesBreadthFirstFinish() const 
{
	debug_func_cerr( 5 );
	
  assert( 0 );
  return NULL;
}

//---------------------------------------------------------< Input/Output routines >---------------------------------------------
void ImplTree::Write( std::ostream& output ) const 
{
	debug_func_cerr( 5 );
	
  for (Node n =0; n < mCurrentNode; n++) {
    std::cout << n << " " << mTree[n] << endl;
  }

}         




} /* namespace alignlib */

