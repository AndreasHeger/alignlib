/*
  alignlib - a library for aligning protein sequences

  $Id$

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


#ifndef ALIGNLIB_FWD_H_
#define ALIGNLIB_FWD_H_

#include "alignlib_types.h"
#include <vector>
#include <map>
#include <string>
#include <boost/shared_ptr.hpp>

namespace alignlib
{
	/** actor objects and their handles */
	class Encoder;
	typedef boost::shared_ptr<Encoder>HEncoder;
	
	class Weightor;
	typedef boost::shared_ptr<Weightor>HWeightor;
	
	class Regularizor;
	typedef boost::shared_ptr<Regularizor>HRegularizor;
	
	class LogOddor;
	typedef boost::shared_ptr<LogOddor>HLogOddor;
	
	class MultipleAlignment;
	typedef boost::shared_ptr<MultipleAlignment>HMultipleAlignment;
	
	class Renderer;
	typedef boost::shared_ptr<Renderer>HRenderer;	
	
	class Alignment;
	typedef boost::shared_ptr<Alignment>HAlignment;
	
	class AlignmentIterator;
	typedef boost::shared_ptr<AlignmentIterator>HAlignmentIterator;
	
	class Alignatum;
	typedef boost::shared_ptr<Alignatum>HAlignatum;
	
	class Alignandum;
	typedef boost::shared_ptr<Alignandum>HAlignandum;

	class Sequence;
	typedef boost::shared_ptr<Sequence>HSequence;

	class Profile;
	typedef boost::shared_ptr<Profile>HProfile;

	class Alignator;
	typedef boost::shared_ptr<Alignator>HAlignator;
	
	class Scorer;
	typedef boost::shared_ptr<Scorer>HScorer;
	
	class Treetor;
	typedef boost::shared_ptr<Treetor>HTreetor;
	
	class Tree;
	typedef boost::shared_ptr<Tree>HTree;
	
	class Distor;
	typedef boost::shared_ptr<Distor>HDistor;

	class DistanceMatrix;
	typedef boost::shared_ptr<DistanceMatrix>HDistanceMatrix;
	
	class Fragmentor;
	typedef boost::shared_ptr<Fragmentor>HFragmentor;
	
	class Iterator2D;
	typedef boost::shared_ptr<Iterator2D>HIterator2D;
	
	/** various matrix definitions */
	template<class T> class Matrix;
	
	typedef Matrix<Score> ScoreMatrix;
    typedef boost::shared_ptr<ScoreMatrix>HScoreMatrix;
	  
	typedef Matrix<Frequency> FrequencyMatrix;
    typedef boost::shared_ptr<FrequencyMatrix>HFrequencyMatrix;
	
	typedef Matrix<Count> CountMatrix;
    typedef boost::shared_ptr<CountMatrix>HCountMatrix;
	
    typedef Matrix<Score> MutationMatrix;
    typedef boost::shared_ptr<MutationMatrix>HMutationMatrix;
    
    typedef Matrix<Score> SubstitutionMatrix;	
    typedef boost::shared_ptr<SubstitutionMatrix>HSubstitutionMatrix;
    
    /** A vector of Residues */
    typedef std::vector< Residue > ResidueVector;
    typedef boost::shared_ptr<ResidueVector>HResidueVector;

    /** A vector of Frequencies */
    typedef std::vector< Frequency> FrequencyVector;
    typedef boost::shared_ptr<FrequencyVector>HFrequencyVector;    
    
    /** A vector of aligned fragmens */
    typedef std::vector<HAlignment> FragmentVector;
    typedef boost::shared_ptr<FragmentVector>HFragmentVector;    

    /** A vector of positions */
    typedef std::vector< Position > PositionVector;
    typedef boost::shared_ptr<PositionVector>HPositionVector;    

    /** A vector of scores */
    typedef std::vector< Score > ScoreVector;
    typedef boost::shared_ptr<ScoreVector>HScoreVector;    

    /** A vector with sequence weights
     */
    typedef std::vector< SequenceWeight > SequenceWeights;
    typedef boost::shared_ptr<SequenceWeights>HSequenceWeights;    

    /** A vector with nodes
     */
    typedef std::vector<Node> NodeVector;
    typedef boost::shared_ptr<NodeVector>HNodeVector;
    
    /** A palette for colouring residues */
    typedef std::map<unsigned char,std::string>Palette;
    typedef boost::shared_ptr<Palette>HPalette;    

}



#endif /*ALIGNLIB_FWD_H_*/
