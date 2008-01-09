/*
  alignlib - a library for aligning protein sequences

  $Id: test_Alignator.cpp,v 1.3 2004/06/02 12:11:38 aheger Exp $

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

/** Test alignata objects

    note: for AlignatorIdentity, AlignatorSimilarity, etc. the 
    alignments look terrible, since they are wraparound pairwise
    alignments.

 */

#if HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>
#include <fstream>
#include <vector>


#include <time.h> 
#include <assert.h>

#include "alignlib.h"

#include "Alignandum.h"
#include "HelpersSequence.h"

#include "Iterator2D.h"
#include "HelpersIterator2D.h"

using namespace std;
using namespace alignlib;

void Print( HIterator2D & iterator )
{
    {
      Iterator2D::const_iterator ir( iterator->row_begin()), ir_end( iterator->row_end() );
      Iterator2D::const_iterator ic( iterator->col_begin()), ic_end( iterator->col_end() );
      std::cout << "iterator range: row=" << *ir << " " << *ir_end << " col=" << *ic << " " << *ic_end << std::endl;
    }

        
    Iterator2D::const_iterator ir( iterator->row_begin()), ir_end( iterator->row_end() );
	for (; ir != ir_end; ++ir)
	{
		std::cout << *ir << " : ";
		Iterator2D::const_iterator ic( iterator->col_begin( *ir )), ic_end( iterator->col_end( *ir ) );
		for (; ic != ic_end; ++ic)
		{
			std::cout << "\t" << *ic;
		}
		std::cout << std::endl;
	}

	Iterator2D::const_iterator ic( iterator->col_begin()), ic_end( iterator->col_end() );
	for (; ic != ic_end; ++ic)
	{
		std::cout << *ic << " : ";
		Iterator2D::const_iterator ir( iterator->row_begin( *ic )), ir_end( iterator->row_end( *ic ) );
		for (; ir != ir_end; ++ir)
		{
			std::cout << "\t" << *ir;
		}
		std::cout << std::endl;
	}

	
}


int main () {

	HAlignandum seq1(makeSequence( "ACDEFGHIKL"));
	HAlignandum seq2(makeSequence( "ACDEFGHIKLMN"));	  

	{
		std::cout << "--------------------- testing Iterator2DFull ----------------------------------" << std::endl;
		HIterator2D iterator(makeIterator2DFull( seq1, seq2));
		assert( iterator->row_size() == seq1->getLength() ); 
		assert( iterator->col_size() == seq2->getLength() );
		Print( iterator );        
	}

	{
		std::cout << "--------------------- testing Iterator2DFull with useSegment ----------------------------------" << std::endl;
		seq1->useSegment( 3, 6);
		HIterator2D iterator(makeIterator2DFull( seq1, seq2));
		assert( iterator->row_size() == 3 );
                assert( iterator->row_front() == 3 );
                assert( iterator->row_back() == 5 );                
		assert( iterator->col_size() == seq2->getLength() );
		Print( iterator );        
		seq1->useSegment();
		seq2->useSegment( 3, 6);
		iterator = makeIterator2DFull( seq1, seq2);
		assert( iterator->row_size() == seq1->getLength() ); 
		assert( iterator->col_size() == 3);
                assert( iterator->col_front() == 3 );
                assert( iterator->col_back() == 5 );                		
		Print( iterator );        
	
		  {
		    HIterator2D iterator2(iterator->getNew( seq1, seq2));
	            assert( iterator2->row_size() == seq1->getLength() ); 
	                assert( iterator2->col_size() == 3);
	                assert( iterator2->col_front() == 3 );
	                assert( iterator2->col_back() == 5 );
		  }
		
		  seq2->useSegment();
	}


	{
		std::cout << "--------------------- testing Iterator2DBanded with diagonals -2, 1 -----" << std::endl;    
		HIterator2D iterator = makeIterator2DBanded( seq1, seq2, -2, 1);
		std::cout << "lengths" << std::endl;
		std::cout << iterator->row_size() << iterator->row_front() << " " << iterator->row_back() << std::endl;
		std::cout << seq1->getLength() << std::endl;
		assert( iterator->row_size() == seq1->getLength() ); 
		assert( iterator->col_size() == seq2->getLength() -1);
		Print( iterator );
	}

	{
		std::cout << "--------------------- testing Iterator2DBanded with diagonals -4, -2 -----" << std::endl;    
		HIterator2D iterator(makeIterator2DBanded( seq1, seq2, -4, -2));
		assert( iterator->row_size() == seq1->getLength() -2); 
		assert( iterator->col_size() == seq2->getLength() -4 );
		Print( iterator );
	}

	{
		std::cout << "--------------------- testing Iterator2DBanded with diagonals 2, 4 -----" << std::endl;    
		HIterator2D iterator = makeIterator2DBanded( seq1, seq2, 2, 4);
		assert( iterator->row_size() == seq1->getLength() ); 
		assert( iterator->col_size() == seq2->getLength() -2);
		Print( iterator );
	}
}
