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
#include "HelpersSequence.h"
#include "HelpersProfile.h"
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

	MultipleAlignment * reference = makeMultipleAlignment();

	{
		reference->add(makeAlignatumFromString("0123456789"));
		reference->add(makeAlignatumFromString("0123456789"));
		reference->add(makeAlignatumFromString("0123456789"));
	}

	std::auto_ptr<Alignment> ali(makeAlignmentSet());
	{
		ali->addPair( new ResiduePAIR( 2,2, 1));
		ali->addPair( new ResiduePAIR( 3,3, 1));
		ali->addPair( new ResiduePAIR( 5,4, 1));
		ali->addPair( new ResiduePAIR( 6,5, 1));
		ali->addPair( new ResiduePAIR( 7,7, 1));
		ali->addPair( new ResiduePAIR( 8,8, 1));

		cout << *ali << endl;

	}	

	{
		// create a multiple alignment 
		std::auto_ptr<MultipleAlignment>m1(makeMultipleAlignment());

		m1->add(makeAlignatumFromString("-AAAAA-CCAAA-"));
		m1->add(makeAlignatumFromString("AAAAAA-CCAAAA"));
		m1->add(makeAlignatumFromString("-A-AAA-CCA-A-"));
		m1->add(makeAlignatumFromString("AAAAAAA--AAAA"));

		assert( m1->getLength() == 13);
		assert( m1->getWidth() == 4);

		std::cout << "testing ouptut" << std::endl;
		cout << *m1 << endl;
		// write segments of the alignment and check numbering
		// no change
		m1->write( cout, 0, 13 ); cout << endl;
		// truncated
		m1->write( cout, 1, 12 ); cout << endl;
		// truncated
		m1->write( cout, 2, 11 ); cout << endl;
		//truncated
		m1->write( cout, 3, 10 ); cout << endl;
	}

	{
		std::auto_ptr<MultipleAlignment> m1(reference->getClone());
		std::cout << "testing using sstream" << std::endl;
		std::ostringstream temp;
		std::ostream * p = &temp;
		m1->write( *p ); *p << ends; std::cout << temp.str() <<endl;

	}

	{
		std::cout << "testing directly adding a multiple alignment" << std::endl;
		std::auto_ptr<MultipleAlignment> m1(reference->getClone());
		std::auto_ptr<MultipleAlignment> m2(makeMultipleAlignment());

		m2->add( &*m1 );
		m2->add( &*m1 );

		assert( m2->getWidth() == m1->getWidth() * 2);
		assert( m2->getLength() == m2->getLength() );
		cout << *m2 << endl;
	}

	{
		std::cout << "testing adding alignments with alignment between them" << std::endl;
		std::auto_ptr<MultipleAlignment> m1(reference->getClone());
		std::auto_ptr<MultipleAlignment> m2(m1->getClone());

		m2->add( &*m1, &*ali, true );
		assert( m2->getLength() == 8);
		assert( m2->getWidth() == 2 * m1->getWidth() );

		debug_cerr(5, "mali is in row\n" << *m2 );
	}

	{
		std::auto_ptr<MultipleAlignment> m1(reference->getClone());
		std::auto_ptr<MultipleAlignment> m2(m1->getClone());

		m2->add( &*m1, &*ali, false);

		assert( m2->getLength() == 8);
		assert( m2->getWidth() == 2 * m1->getWidth() );

		debug_cerr(5, "mali is in col\n" << *m2 );
	}

	//-------------------------add new objects using an aligment 
	{	
		std::auto_ptr<MultipleAlignment> m1(reference->getClone());

		m1->add(makeAlignatumFromString("0123456789"), &*ali, true );
		cout << *m1 << endl;
	}	

	{
		std::auto_ptr<MultipleAlignment> m1(reference->getClone());

		m1->add(makeAlignatumFromString("0123456789"), &*ali, false );
		cout << *m1 << endl;
	}	
		
	// check mali dots
	std::cout << "## checking MultipleAlignmentDots" << std::endl;
	{

		// create a multiple alignment 
		MultipleAlignment * m1 = makeMultipleAlignmentDots( false );
		MultipleAlignment * m2 = makeMultipleAlignmentDots( false );    

		Alignment * a1 = makeAlignmentVector();
		fillAlignmentIdentity( a1, 1, 4, 0 );
		fillAlignmentIdentity( a1, 6, 8, 0 );

		Alignment * a2 = makeAlignmentVector();
		fillAlignmentIdentity( a2, 0, 5, +1 );

		Alignment * a3 = makeAlignmentVector();
		fillAlignmentIdentity( a3, 1, 6, -1 );
		fillAlignmentIdentity( a3, 6, 7, 0 );
		fillAlignmentIdentity( a3, 7, 8, +1 );
		
		m1->add(makeAlignatumFromString("ABCDGHJL"), a1);
		m1->add(makeAlignatumFromString(".ABCDEFGH"), a2);
		m1->add(makeAlignatumFromString("BCDEFIJKLM"), a3);

		m2->add(makeAlignatumFromString("ABCDGHJL"), a1);
		m2->add(makeAlignatumFromString(".ABCDEFGH"), a2);
		m2->add(makeAlignatumFromString("BCDEFIJKLM"), a3);

		cout << *m1 << endl;
		cout << *m2 << endl;    

		Alignment * aa1 = makeAlignmentVector();
		fillAlignmentIdentity( aa1, 1, 9, 0 );    

		// IMPORTANT: m1 takes ownership of objects in m2.
		// Thus: do not delete m2.
		m1->add( m2, aa1 );

		cout << *m1 << endl;        
		delete a1;
		delete a2;
		delete a3;
		delete m1;
		delete m2;
		delete aa1;

	}

	{

		// create a multiple alignment 
		MultipleAlignment * m1 = makeMultipleAlignment();
		MultipleAlignment * m2 = makeMultipleAlignmentDots( false );    

		Alignment * a1 = makeAlignmentVector();
		fillAlignmentIdentity( a1, 0, 8, 0 );

		m1->add(makeAlignatumFromString("ABCDGHIJL"), a1);
		m1->add(makeAlignatumFromString("ABCDGHIJL"), a1);

		m2->add(makeAlignatumFromString("ABCDGHIJL"), a1);
		m2->add(makeAlignatumFromString("ABCDGHIJL"), a1);    

		cout << *m1 << endl;
		cout << *m2 << endl;    

		Alignment * aa1 = makeAlignmentVector();
		fillAlignmentIdentity( aa1, 0, 5 );
		fillAlignmentIdentity( aa1, 5, 7, 1 );        

		m1->add( m2, aa1 );

		cout << *m1 << endl;

		delete m1;
		delete m2;
		delete a1;
		delete aa1;
	}

	// check mali dots
	std::cout << "## checking MultipleAlignmentDots options" << std::endl;
	{

		// create a multiple alignment 
		MultipleAlignment * m1 = makeMultipleAlignmentDots( true );

		Alignment * a1 = makeAlignmentVector();

		Alignment * a2 = makeAlignmentVector();
		fillAlignmentIdentity( a2, 0, 5, +1 );

		Alignment * a3 = makeAlignmentVector();
		fillAlignmentIdentity( a3, 1, 6, -1 );
		fillAlignmentIdentity( a3, 6, 7, 0 );
		fillAlignmentIdentity( a3, 7, 8, +1 );

		m1->add(makeAlignatumFromString("ABCDGHJL") );
		m1->add(makeAlignatumFromString("YABCDEFGH"), a2);
		m1->add(makeAlignatumFromString("BCDEFIJKLM"), a3);

		std::cout << *m1 << endl;

		delete m1;
		delete a1;
		delete a2;
		delete a3;

	}

	{
		std::cout << "## checking MultipleAlignmentDots options" << std::endl;

		// create a multiple alignment 
		MultipleAlignment * m1 = makeMultipleAlignmentDots( true );

		Alignment * a1 = makeAlignmentVector();
		fillAlignmentIdentity( a1, 1, 4, 0 );
		fillAlignmentIdentity( a1, 6, 8, 0 );

		Alignment * a2 = makeAlignmentVector();
		fillAlignmentIdentity( a2, 0, 5, +1 );

		Alignment * a3 = makeAlignmentVector();
		fillAlignmentIdentity( a3, 1, 6, -1 );
		fillAlignmentIdentity( a3, 6, 7, 0 );
		fillAlignmentIdentity( a3, 7, 8, +1 );

		m1->add(makeAlignatumFromString("ABCDGHJL"), a1);
		m1->add(makeAlignatumFromString("YABCDEFGH"), a2);
		m1->add(makeAlignatumFromString("BCDEFIJKLM"), a3);

		std::cout << *m1 << endl;

		delete m1;
		delete a1;
		delete a2;
		delete a3;

	}

	// check mali dots

	{
		std::cout << "## checking MultipleAlignmentDots options" << std::endl;

		// create a multiple alignment 
		MultipleAlignment * m1 = makeMultipleAlignmentDots( true, 0 );

		Alignment * a1 = makeAlignmentVector();
		fillAlignmentIdentity( a1, 1, 4, 0 );
		fillAlignmentIdentity( a1, 6, 8, 0 );

		Alignment * a2 = makeAlignmentVector();
		fillAlignmentIdentity( a2, 0, 5, +1 );

		Alignment * a3 = makeAlignmentVector();
		fillAlignmentIdentity( a3, 1, 6, -1 );
		fillAlignmentIdentity( a3, 6, 7, 0 );
		fillAlignmentIdentity( a3, 7, 8, +1 );

		m1->add(makeAlignatumFromString("ABCDGHJL"), a1);
		m1->add(makeAlignatumFromString(".ABCDEFGH"), a2);
		m1->add(makeAlignatumFromString("BCDEFIJKLM"), a3);

		cout << *m1 << endl;

		delete m1;
		delete a1;
		delete a2;
		delete a3;
	}
	
	delete reference;
}







