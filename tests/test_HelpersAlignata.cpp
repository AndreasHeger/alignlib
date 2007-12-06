/*
  alignlib - a library for aligning protein sequences

  $Id: test_HelpersAlignata.cpp,v 1.5 2004/10/14 23:34:09 aheger Exp $

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
 */

#if HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>
#include <fstream>
#include <sstream>

#include <time.h> 

#include "alignlib.h"
#include "Alignandum.h"
#include "HelpersSequence.h"

#include "Alignata.h"
#include "HelpersAlignata.h"
#include "AlignataIterator.h"

using namespace std;

#ifdef WITH_DMALLOC
#include <dmalloc.h>
#endif

using namespace alignlib;

bool isIdentical( const Alignata * a, const Alignata * b, bool inverse = false) {

	AlignataConstIterator it1(a->begin());
	AlignataConstIterator it1_end(a->end());

	AlignataConstIterator it2(b->begin());
	AlignataConstIterator it2_end(b->end());

	bool is_identical = true;

	for (; it1 != it1_end; ++it1, ++it2) {
		if (!inverse) {
			if (it1->mRow != it2->mRow && it1->mCol != it2->mCol) 
				is_identical = false;
		} else {
			if (it1->mRow != it2->mCol && it1->mCol != it2->mRow) 
				is_identical = false;
		}
	}

	return is_identical;
}

bool TestCompressionMonotone( Alignata * a, 
		const char *xrow, 
		const char *xcol,
		const char *crow = NULL, 
		const char *ccol = NULL) {

	if (crow == NULL) crow = xrow;
	if (ccol == NULL) ccol = xcol;

	std::string row (xrow);
	std::string col (xcol);

	a = fillAlignataCompressed( a, 3, row, 3, col);

	std::stringstream output;

	writeAlignataCompressed( output, a  );
	output.seekp( 0);
	std::string new_row;
	std::string new_col;
	
	output >> new_row >> new_col;
	
	row = crow;
	col = ccol;

	if (row == new_row && col == new_col) 
		return true;
	else {
		cout << "xrow=" << xrow << " xcol=" << xcol << endl;
		cout << "crow=" << crow << " ccol=" << ccol << endl;
		cout << *a << endl;
		cout << "row=" << row << " col=" << col << " new_row=" << new_row << " new_col=" << new_col << endl;
		return false;
	}
}

bool TestCompressionDiagonal( Alignata * a, 
		const char *xrow, 
		const char *crow = NULL) {


	if (crow == NULL) crow = xrow;

	std::string row (xrow);

	a = fillAlignataCompressedDiagonal( a, row );

	std::stringstream output;
	writeAlignataCompressedDiagonal( output, a );
	std::string new_row(output.str());

	row = crow;

	if (row == new_row) 
		return true;
	else {
		cout << "xrow=" << xrow << " crow=" << crow << endl;
		cout << *a << endl;
		cout << "new_row=" << new_row << endl;
		return false;
	}
}


// tests for diagonal alignments
void TestDiagonal(Alignata * a) 
{

	{ 
		cout << "testing...fill/writeAlignataCompressedDiagonal()...";

		if (TestCompressionDiagonal( a, "-1:-5+2" ) &&				
				TestCompressionDiagonal( a, "-3:-5+2-2+2;5:-2+3" ) &&				
				TestCompressionDiagonal( a, "3:-5+2" ) &&
				TestCompressionDiagonal( a, "-3:-5+2;0:-5+4" )
		) 
			cout << "passed" << endl;
		else
			cout << "failed" << endl;      
	}

}


// tests for monotone alignments
void TestMonotone( Alignata * a) {

	{ 
		cout << "testing...fill/writeAlignataCompressed()...";

		if (TestCompressionMonotone( a, "+5-5+5", "+15") &&				// treat gaps in one sequence correctly
				TestCompressionMonotone( a, "+15", "+5-5+5") &&				// treat gaps in one sequence correctly
				TestCompressionMonotone( a, "+5-5+5", "+5-5+5", "+10", "+10") &&	// skip a gap in both sequences
				TestCompressionMonotone( a, "-5+10", "+10-5", "+5", "+5") &&		// skip end gaps
				TestCompressionMonotone( a, "+10-5", "-5+10", "+5", "+5") &&		// skip end gaps
				TestCompressionMonotone( a, "-5+5-5+10", "-5+5-5+5-5", "+10", "+10") )	// skip end and middle gaps gaps
			cout << "passed" << endl;
		else
			cout << "failed" << endl;      
	}

	{
		cout << "testing fillAlignataSummation()...";
		a->clear();
		a->addPair( new ResiduePAIR( 3, 4, 0));
		a->addPair( new ResiduePAIR( 4, 5, 0));         
		a->addPair( new ResiduePAIR( 5, 7, 0));                  
		a->addPair( new ResiduePAIR( 9, 9, 0));                           

		Alignata * b = a->getNew();
		Alignata * c = a->getNew();

		fillAlignataSummation( b, c, a, true, true, true, true, 11, 11);
	}
}

// tests for both
void Test( Alignata * a) {
	{ 
		cout << "testing...copyAlignata()...";

		Alignata * a_new = makeAlignataVector();
		copyAlignata( a_new, a);
		delete a_new;
	}

	{ 
		cout << "testing...combineAlignata()...";
		Alignata * a_new = makeAlignataVector();
		combineAlignata( a_new, a, a, alignlib::RR); 
		combineAlignata( a_new, a, a, alignlib::CR);
		combineAlignata( a_new, a, a, alignlib::RC);
		combineAlignata( a_new, a, a, alignlib::CC);
		delete a_new;

		cout << "passed" << endl;
	}

	{
		cout << "testing...rescoreAlignment()...";
		rescoreAlignment( a, 1 );
		cout << "passed" << endl;
	}
	{
		cout << "testing...rescoreAlignment()...";
		Alignandum * s1 = makeSequence("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA");
		Alignandum * s2 = makeSequence("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA");
		rescoreAlignmentPrivate( a, s1, s2);
		delete s1;
		delete s2;
		cout << "passed" << endl;
	}

	{
		cout << "testing...fillAlignataIdentity()...";
		fillAlignataIdentity( a, 1, 10, 0);
		cout << "passed" << endl;
	}

	{ 
		cout << "testing...writeAlignataTable()...";
		ostringstream out;
		writeAlignataTable( out, a);
		cout << "passed" << endl;
	}

	{ 
		cout << "testing...addAlignata2Alignata()...";
		Alignata * b = makeAlignataVector();
		addAlignata2Alignata( b, a );
		delete b;
		cout << "passed" << endl;
	}

	{ 
		cout << "testing...addMappedAlignata2Alignata()...";
		Alignata * b = makeAlignataVector();
		addMappedAlignata2Alignata( b, a, a, alignlib::RR );
		delete b;
		cout << "passed" << endl;
	}

}



int main () {

	Alignata * a;

	cout << "---------------------Testing AlignataVector-------------------------------" << endl;
	a = makeAlignataVector();
	TestMonotone( a );
	Test( a );
	delete a;

	cout << "---------------------Testing AlignataSet----------------------------------" << endl;
	a = makeAlignataSet();
	TestMonotone( a );
	Test( a );
	delete a;

	cout << "---------------------Testing AlignataHash----------------------------------" << endl;
	a = makeAlignataHash();
	TestMonotone( a );
	Test( a );
	delete a;

	cout << "---------------------Testing AlignataSetCol------------------------------" << endl;
	a = makeAlignataSetCol();
	TestMonotone( a );
	Test( a );
	delete a;

	cout << "---------------------Testing AlignataMatrixRow-------------------------------" << endl;
	a = makeAlignataMatrixRow();
	TestMonotone( a );
	Test( a );
	delete a;

	cout << "---------------------Testing AlignataMatrixDiagonal-------------------------------" << endl;
	a = makeAlignataMatrixDiagonal();
	TestDiagonal( a );
	Test( a );
	delete a;

	cout << "---------------------Testing AlignataHashDiagonal------------------------------" << endl;
	a = makeAlignataHashDiagonal();
	TestDiagonal( a );
	Test( a );
	delete a;

	cout << "---------------------Testing AlignataMatrixUnsorted-------------------------------" << endl;
	a = makeAlignataMatrixUnsorted();
	TestDiagonal( a );
	Test( a );
	delete a;

}








