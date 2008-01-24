/*
  alignlib - a library for aligning protein sequences

  $Id: test_MultipleAlignment.cpp,v 1.6 2004/06/02 12:11:38 aheger Exp $

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

/** Test the MultipleAlignment - object
 */

#if HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>
#include <fstream>
#include <sstream>
#include <cassert>

#include "Alignandum.h"
#include "HelpersAlignandum.h"
#include "Alignator.h"
#include "HelpersAlignator.h"
#include "Alignatum.h"
#include "HelpersAlignatum.h"
#include "Alignment.h"
#include "HelpersAlignment.h"
#include "MultipleAlignment.h"
#include "HelpersMultipleAlignment.h"
#include "AlignlibDebug.h"

using namespace std;
using namespace alignlib;

const char FILE_SEQ[] = "data/test.seq";
const char FILE_FASTA[] = "data/test.fasta";


int main () 
{

	HMultipleAlignment reference(makeMultipleAlignment());

	{
		reference->add(makeAlignatum("0123456789"));
		reference->add(makeAlignatum("0123456789"));
		reference->add(makeAlignatum("0123456789"));
	}

	HAlignment ali(makeAlignmentSet());
	{
		ali->addPair( ResiduePair( 2,2, 1));
		ali->addPair( ResiduePair( 3,3, 1));
		ali->addPair( ResiduePair( 5,4, 1));
		ali->addPair( ResiduePair( 6,5, 1));
		ali->addPair( ResiduePair( 7,7, 1));
		ali->addPair( ResiduePair( 8,8, 1));

		cout << *ali << endl;

	}	

	{
		// create a multiple alignment 
		HMultipleAlignment m1(makeMultipleAlignment());

		m1->add(makeAlignatum("-AAAAA-CCAAA-"));
		m1->add(makeAlignatum("AAAAAA-CCAAAA"));
		m1->add(makeAlignatum("-A-AAA-CCA-A-"));
		m1->add(makeAlignatum("AAAAAAA--AAAA"));

		assert( m1->getLength() == 13);
		assert( m1->getNumSequences() == 4);

		std::cout << "testing ouptut" << std::endl;
		cout << *m1 << endl;
	}

	{
		HMultipleAlignment m1(reference->getClone());
		std::cout << "testing using sstream" << std::endl;
		std::ostringstream temp;
		std::ostream * p = &temp;
		m1->write( *p ); *p << ends; std::cout << temp.str() <<endl;

	}

	{
		std::cout << "testing directly adding a multiple alignment" << std::endl;
		HMultipleAlignment m1(reference->getClone());
		HMultipleAlignment m2(makeMultipleAlignment());

		m2->add( m1 );
		m2->add( m1 );

		assert( m2->getNumSequences() == m1->getNumSequences() * 2);
		assert( m2->getLength() == m2->getLength() );
		cout << *m2 << endl;
	}

	{
		std::cout << "testing adding alignments with alignment between them" << std::endl;
		HMultipleAlignment m1(reference->getClone());
		HMultipleAlignment m2(m1->getClone());

		m2->add( m1, ali, true );
		assert( m2->getLength() == 8);
		assert( m2->getNumSequences() == 2 * m1->getNumSequences() );

		debug_cerr(5, "mali is in row\n" << *m2 );
	}

	{
		HMultipleAlignment m1(reference->getClone());
		HMultipleAlignment m2(m1->getClone());

		m2->add( m1, ali, false);

		assert( m2->getLength() == 8);
		assert( m2->getNumSequences() == 2 * m1->getNumSequences() );

		debug_cerr(5, "mali is in col\n" << *m2 );
	}

	//-------------------------add new objects using an aligment 
	{	
		HMultipleAlignment m1(reference->getClone());

		m1->add(makeAlignatum("0123456789"), ali, true );
		cout << *m1 << endl;
	}	

	{
		HMultipleAlignment m1(reference->getClone());

		m1->add(makeAlignatum("0123456789"), ali, false );
		cout << *m1 << endl;
	}	
		
	// check mali dots
	std::cout << "## checking MultipleAlignmentDots" << std::endl;
	{

		// create a multiple alignment 
		HMultipleAlignment m1(makeMultipleAlignmentDots( false ));
		HMultipleAlignment m2(makeMultipleAlignmentDots( false ));    

		HAlignment a1(makeAlignmentVector());
		addDiagonal2Alignment( a1, 1, 4, 0 );
		addDiagonal2Alignment( a1, 6, 8, 0 );

		HAlignment a2(makeAlignmentVector());
		addDiagonal2Alignment( a2, 0, 5, +1 );

		HAlignment a3(makeAlignmentVector());
		addDiagonal2Alignment( a3, 1, 6, -1 );
		addDiagonal2Alignment( a3, 6, 7, 0 );
		addDiagonal2Alignment( a3, 7, 8, +1 );
		
		m1->add(makeAlignatum("ABCDGHJL"), a1);
		m1->add(makeAlignatum(".ABCDEFGH"), a2);
		m1->add(makeAlignatum("BCDEFIJKLM"), a3);

		m2->add(makeAlignatum("ABCDGHJL"), a1);
		m2->add(makeAlignatum(".ABCDEFGH"), a2);
		m2->add(makeAlignatum("BCDEFIJKLM"), a3);

		cout << *m1 << endl;
		cout << *m2 << endl;    

		HAlignment aa1(makeAlignmentVector());
		addDiagonal2Alignment( aa1, 1, 9, 0 );    
		m1->add( m2, aa1 );
		cout << *m1 << endl;        
	}

	{

		// create a multiple alignment 
		HMultipleAlignment m1(makeMultipleAlignment());
		HMultipleAlignment m2(makeMultipleAlignmentDots( false ));    

		HAlignment a1(makeAlignmentVector());
		addDiagonal2Alignment( a1, 0, 8, 0 );

		m1->add(makeAlignatum("ABCDGHIJL"), a1);
		m1->add(makeAlignatum("ABCDGHIJL"), a1);

		m2->add(makeAlignatum("ABCDGHIJL"), a1);
		m2->add(makeAlignatum("ABCDGHIJL"), a1);    

		cout << *m1 << endl;
		cout << *m2 << endl;    

		HAlignment aa1(makeAlignmentVector());
		addDiagonal2Alignment( aa1, 0, 5 );
		addDiagonal2Alignment( aa1, 5, 7, 1 );        

		m1->add( m2, aa1 );

		cout << *m1 << endl;

	}

	// check mali dots
	std::cout << "## checking MultipleAlignmentDots options" << std::endl;
	{

		// create a multiple alignment 
		HMultipleAlignment m1(makeMultipleAlignmentDots( true ));

		HAlignment a1(makeAlignmentVector());

		HAlignment a2(makeAlignmentVector());
		addDiagonal2Alignment( a2, 0, 5, +1 );

		HAlignment a3(makeAlignmentVector());
		addDiagonal2Alignment( a3, 1, 6, -1 );
		addDiagonal2Alignment( a3, 6, 7, 0 );
		addDiagonal2Alignment( a3, 7, 8, +1 );

		m1->add(makeAlignatum("ABCDGHJL") );
		m1->add(makeAlignatum("YABCDEFGH"), a2);
		m1->add(makeAlignatum("BCDEFIJKLM"), a3);

		std::cout << *m1 << endl;

	}

	{
		std::cout << "## checking MultipleAlignmentDots options" << std::endl;

		// create a multiple alignment 
		HMultipleAlignment m1(makeMultipleAlignmentDots( true ));

		HAlignment a1(makeAlignmentVector());
		addDiagonal2Alignment( a1, 1, 4, 0 );
		addDiagonal2Alignment( a1, 6, 8, 0 );

		HAlignment a2(makeAlignmentVector());
		addDiagonal2Alignment( a2, 0, 5, +1 );

		HAlignment a3(makeAlignmentVector());
		addDiagonal2Alignment( a3, 1, 6, -1 );
		addDiagonal2Alignment( a3, 6, 7, 0 );
		addDiagonal2Alignment( a3, 7, 8, +1 );

		m1->add(makeAlignatum("ABCDGHJL"), a1);
		m1->add(makeAlignatum("YABCDEFGH"), a2);
		m1->add(makeAlignatum("BCDEFIJKLM"), a3);

		std::cout << *m1 << endl;
	}

	// check mali dots

	{
		std::cout << "## checking MultipleAlignmentDots options" << std::endl;

		// create a multiple alignment 
		HMultipleAlignment m1(makeMultipleAlignmentDots( true, 0 ));

		HAlignment a1(makeAlignmentVector());
		addDiagonal2Alignment( a1, 1, 4, 0 );
		addDiagonal2Alignment( a1, 6, 8, 0 );

		HAlignment a2(makeAlignmentVector());
		addDiagonal2Alignment( a2, 0, 5, +1 );

		HAlignment a3(makeAlignmentVector());
		addDiagonal2Alignment( a3, 1, 6, -1 );
		addDiagonal2Alignment( a3, 6, 7, 0 );
		addDiagonal2Alignment( a3, 7, 8, +1 );

		m1->add(makeAlignatum("ABCDGHJL"), a1);
		m1->add(makeAlignatum(".ABCDEFGH"), a2);
		m1->add(makeAlignatum("BCDEFIJKLM"), a3);

		cout << *m1 << endl;

	}
	
}







