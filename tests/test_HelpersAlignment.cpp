/*
  alignlib - a library for aligning protein sequences

  $Id: test_HelpersAlignment.cpp,v 1.5 2004/10/14 23:34:09 aheger Exp $

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

#include "Alignment.h"
#include "HelpersAlignment.h"
#include "AlignmentIterator.h"

using namespace std;

#ifdef WITH_DMALLOC
#include <dmalloc.h>
#endif

using namespace alignlib;

bool isIdentical( const Alignment * a, const Alignment * b, bool inverse = false) {

	AlignmentConstIterator it1(a->begin());
	AlignmentConstIterator it1_end(a->end());

	AlignmentConstIterator it2(b->begin());
	AlignmentConstIterator it2_end(b->end());

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

bool TestCompressionMonotone( Alignment * a, 
		const char *xrow, 
		const char *xcol,
		const char *crow = NULL, 
		const char *ccol = NULL) {

	if (crow == NULL) crow = xrow;
	if (ccol == NULL) ccol = xcol;

	std::string row (xrow);
	std::string col (xcol);

	a = fillAlignmentCompressed( a, 3, row, 3, col);

	std::stringstream output;

	writeAlignmentCompressed( output, a  );
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

bool TestCompressionDiagonal( Alignment * a, 
		const char *xrow, 
		const char *crow = NULL) {


	if (crow == NULL) crow = xrow;

	std::string row (xrow);

	a = fillAlignmentCompressedDiagonal( a, row );

	std::stringstream output;
	writeAlignmentCompressedDiagonal( output, a );
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
void TestDiagonal(Alignment * a) 
{

	{ 
		cout << "testing...fill/writeAlignmentCompressedDiagonal()...";

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
void TestMonotone( Alignment * a) {

	{ 
		cout << "testing...fill/writeAlignmentCompressed()...";

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
		cout << "testing fillAlignmentSummation()...";
		a->clear();
		a->addPair( new ResiduePAIR( 3, 4, 0));
		a->addPair( new ResiduePAIR( 4, 5, 0));         
		a->addPair( new ResiduePAIR( 5, 7, 0));                  
		a->addPair( new ResiduePAIR( 9, 9, 0));                           

		Alignment * b = a->getNew();
		Alignment * c = a->getNew();

		fillAlignmentSummation( b, c, a, true, true, true, true, 11, 11);
	}
}

// tests for both
void Test( Alignment * a) {
	{ 
		cout << "testing...copyAlignment()...";

		Alignment * a_new = makeAlignmentVector();
		copyAlignment( a_new, a);
		delete a_new;
	}

	{ 
		cout << "testing...combineAlignment()...";
		Alignment * a_new = makeAlignmentVector();
		combineAlignment( a_new, a, a, alignlib::RR); 
		combineAlignment( a_new, a, a, alignlib::CR);
		combineAlignment( a_new, a, a, alignlib::RC);
		combineAlignment( a_new, a, a, alignlib::CC);
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
		rescoreAlignment( a, s1, s2);
		delete s1;
		delete s2;
		cout << "passed" << endl;
	}

	{
		cout << "testing...fillAlignmentIdentity()...";
		fillAlignmentIdentity( a, 1, 10, 0);
		cout << "passed" << endl;
	}

	{ 
		cout << "testing...writeAlignmentTable()...";
		ostringstream out;
		writeAlignmentTable( out, a);
		cout << "passed" << endl;
	}

	{ 
		cout << "testing...addAlignment2Alignment()...";
		Alignment * b = makeAlignmentVector();
		addAlignment2Alignment( b, a );
		delete b;
		cout << "passed" << endl;
	}

	{ 
		cout << "testing...addMappedAlignment2Alignment()...";
		Alignment * b = makeAlignmentVector();
		addMappedAlignment2Alignment( b, a, a, alignlib::RR );
		delete b;
		cout << "passed" << endl;
	}

}



int main () {

	Alignment * a;

	cout << "---------------------Testing AlignmentVector-------------------------------" << endl;
	a = makeAlignmentVector();
	TestMonotone( a );
	Test( a );
	delete a;

	cout << "---------------------Testing AlignmentSet----------------------------------" << endl;
	a = makeAlignmentSet();
	TestMonotone( a );
	Test( a );
	delete a;

	cout << "---------------------Testing AlignmentHash----------------------------------" << endl;
	a = makeAlignmentHash();
	TestMonotone( a );
	Test( a );
	delete a;

	cout << "---------------------Testing AlignmentSetCol------------------------------" << endl;
	a = makeAlignmentSetCol();
	TestMonotone( a );
	Test( a );
	delete a;

	cout << "---------------------Testing AlignmentMatrixRow-------------------------------" << endl;
	a = makeAlignmentMatrixRow();
	TestMonotone( a );
	Test( a );
	delete a;

	cout << "---------------------Testing AlignmentMatrixDiagonal-------------------------------" << endl;
	a = makeAlignmentMatrixDiagonal();
	TestDiagonal( a );
	Test( a );
	delete a;

	cout << "---------------------Testing AlignmentHashDiagonal------------------------------" << endl;
	a = makeAlignmentHashDiagonal();
	TestDiagonal( a );
	Test( a );
	delete a;

	cout << "---------------------Testing AlignmentMatrixUnsorted-------------------------------" << endl;
	a = makeAlignmentMatrixUnsorted();
	TestDiagonal( a );
	Test( a );
	delete a;

}








