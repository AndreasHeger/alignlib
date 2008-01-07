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

#include "alignlib.h"
#include <vector>
#include <boost/shared_ptr.hpp>

namespace alignlib
{
	/** actor objects and their handles */
	class Translator;
	typedef boost::shared_ptr<Translator>HTranslator;
	
	class Weightor;
	typedef boost::shared_ptr<Weightor>HWeightor;
	
	class Regularizor;
	typedef boost::shared_ptr<Regularizor>HRegularizor;
	
	class LogOddor;
	typedef boost::shared_ptr<LogOddor>HLogOddor;
	
	class MultipleAlignment;
	typedef boost::shared_ptr<MultipleAlignment>HMultipleAlignment;
	
	class Alignment;
	typedef boost::shared_ptr<Alignment>HAlignment;
	
	class Alignatum;
	typedef boost::shared_ptr<Alignatum>HAlignatum;
	
	class Alignandum;
	typedef boost::shared_ptr<Alignandum>HAlignandum;
	
	class Alignator;
	typedef boost::shared_ptr<Alignator>HAlignator;
	
	class Treetor;
	typedef boost::shared_ptr<Treetor>HTreetor;
	
	class Distor;
	typedef boost::shared_ptr<Distor>HDistor;
	
	class Fragmentor;
	typedef boost::shared_ptr<Fragmentor>HFragmentor;
	
	class Iterator2D;
	typedef boost::shared_ptr<Iterator2D>HIterator2D;
	
	/** various matrix definitions */
	template<class T> class Matrix;
	
	typedef Matrix<Score> ScoreMatrix;
	typedef Matrix<Frequency> FrequencyMatrix;
	typedef Matrix<Count> CountMatrix;
    typedef Matrix<Score> MutationMatrix;
    typedef Matrix<Score> SubstitutionMatrix;	
    
    typedef std::vector< Frequency> FrequencyVector;
    
    
}

#endif /*ALIGNLIB_FWD_H_*/
