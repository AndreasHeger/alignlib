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

#include "alignlib.h"

#include "Alignandum.h"
#include "HelpersSequence.h"

#include "Iterator2D.h"
#include "HelpersIterator2D.h"

#ifdef WITH_DMALLOC
#include <dmalloc.h>
#endif

using namespace std;
using namespace alignlib;

void Print( Iterator2D * iterator )
{
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
}


int main () {

  Alignandum * seq1 = makeSequence( "ACDEFGHIKL");
  {
    std::cout << "--------------------- testing Iterator2DFull ----------------------------------" << std::endl;
    Iterator2D * iterator = makeIterator2DFull( seq1, seq1);
    assert( iterator->row_size() == seq1->getLength() ); 
    assert( iterator->col_size() == seq1->getLength() );
    Print( iterator );    
    delete iterator;
  }
  
  {
    std::cout << "--------------------- testing Iterator2DBanded with diagonals -2, 1 -----" << std::endl;    
    Iterator2D * iterator = makeIterator2DBanded( seq1, seq1, -2, 1);
    assert( iterator->row_size() == seq1->getLength() -2 ); 
    assert( iterator->col_size() == seq1->getLength() );
    Print( iterator );
    delete iterator;
  }
  
  delete seq1;
  
}
