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
#include <sstream>
#include <vector>

#include <time.h> 

#include "alignlib.h"

using namespace std;
using namespace alignlib;

// typedef enum MODE { PAIR, MATRIX };
bool testPairwiseAlignment(int test_id,
		HAlignator & a, 
		const HAlignandum & benchmark_row,
		const HAlignandum & benchmark_col,
		Position row_from,
		Position row_to,
		const char * row,
		Position col_from,
		Position col_to,
		const char * col,
		Score score )
{

	std::cout << "======== test " << test_id << " =========" << std::endl;

	std::cout << *benchmark_row << std::endl;
	std::cout << *benchmark_col << std::endl;

	HAlignment result(makeAlignmentVector());

	a->align( result, benchmark_row, benchmark_col );

	std::string r_row(row);
	std::string r_col(col);

	std::cout << AlignmentFormatExplicit( result, benchmark_row, benchmark_col ) << std::endl ;

	std::cout << "result=" << *result << std::endl;
	
	AlignmentFormatEmissions format( result );
	
	if ( result->getScore() == score &&
			row_from == result->getRowFrom() &&
			row_to == result->getRowTo() &&
			col_from == result->getColFrom() &&
			col_to == result->getColTo() &&
			format.mRowAlignment == r_row &&
			format.mColAlignment == r_col )
	{
		std::cout << "test " << test_id << " success" << std::endl;
		return true;
	}
	else
	{
		std::cout << row_from << "-" << row_to << ":" << r_row
		<< (( format.mRowAlignment == r_row && row_from == result->getRowFrom() && row_to == result->getRowTo() ) ? " == " : " != ")
		<< result->getRowFrom() << "-" << result->getRowTo() << ":" << format.mRowAlignment << "\t"
		<< col_from << "-" << col_to << ":" << r_col
		<< (( format.mColAlignment == r_col && col_from == result->getColFrom() && col_to == result->getColTo() ) ? " == " : " != ") 
		<< result->getColFrom() << "-" << result->getColTo() << ":" << format.mColAlignment << "\t"
		<< score << (( score == result->getScore() ) ? " == " : " != ") << result->getScore()
		<< std::endl;
		std::cout << "test " << test_id << " failure" << std::endl;      

		return false;
	}
}

bool testWrappedAlignment(int test_id,
		HAlignator & a, 
		const HAlignandum & benchmark_row,
		const HAlignandum & benchmark_col,
		const char * r_ali,
		Score score )
{

	std::cout << "======== test " << test_id << " =========" << std::endl;
	std::cout << *benchmark_row << std::endl;
	std::cout << *benchmark_col << std::endl;

	HAlignment result(makeAlignmentMatrixDiagonal());

	a->align( result, benchmark_row, benchmark_col );

	std::string r_row(r_ali);

	AlignmentFormatDiagonals format( result );
	
	
	if ( result->getScore() == score &&
			format.mAlignment == r_ali )
	{
		std::cout << "test " << test_id << " success" << std::endl;
		return true;
	}
	else
	{
		std::cout << r_ali
		<< (( format.mAlignment == r_ali ) ? " == " : " != ")
		<< format.mAlignment
		<< std::endl;
		std::cout << "test " << test_id << " failure" << std::endl;      

		return false;
	}
}


/* test pairwise alignments. For the calculation of the score it is assumed, 
   that the identity matrix is used (match = 1, mismatch = -1, gop = -0.2, gep = -0.1).
 */
/*
void testPairAlignator( Alignator * a,
			Alignment * ali)
{

  { 
    std::cout << "testing...Align( seq1, seq2 )#1...";

    Alignandum * a1 = makeSequence( "AAAAACCCCCAAAAA");
    Alignandum * a2 = makeSequence( "AAAAACCCCCAAAAA");

    if (
	TestPairwiseAlignment( a, ali, a1, a2, "+15", "+15", 15) &&
	TestPairwiseAlignment( a, ali, a2, a1, "+15", "+15", 15) 
	) 
      std::cout << "passed" << std::endl;
    else 
      std::cout << "failed" << std::endl;

    delete a1;
    delete a2;
  }

  { 
    std::cout << "testing...Align( seq1, seq2 )#2...";

    Alignandum * a1 = makeSequence( "AAAAACCCCCAAAAA");
    Alignandum * a2 = makeSequence( "CCCCC");

    if (
	TestPairwiseAlignment( a, ali, a1, a2, "+5", "+5", 5) &&
	TestPairwiseAlignment( a, ali, a2, a1, "+5", "+5", 5) 
	) 
      cout << "passed" << endl;
    else 
      cout << "failed" << endl;

    delete a1;
    delete a2;
  }

  { 
    cout << "testing...Align( seq1, seq2 )#3...";

    Alignandum * a1 = makeSequence( "AAAAACCCAACCCAAAAA");
    Alignandum * a2 = makeSequence( "CCCCCC");

    if (
	TestPairwiseAlignment( a, ali, a1, a2, "+8", "+3-2+3", 4.8) &&
	TestPairwiseAlignment( a, ali, a2, a1, "+3-2+3", "+8", 4.8) 
	) 
      cout << "passed" << endl;
    else 
      cout << "failed" << endl;

    delete a1;
    delete a2;
  }

  { 
    cout << "testing...Align( seq1, seq2 )#4...";

    Alignandum * a1 = makeSequence( "AAAAACCCACCCAAAAA");
    Alignandum * a2 = makeSequence( "CCCKCCC");

    if (
	TestPairwiseAlignment( a, ali, a1, a2, "+7", "+7", 5) &&
	TestPairwiseAlignment( a, ali, a2, a1, "+7", "+7", 5) 
	) 
      cout << "passed" << endl;
    else 
      cout << "failed" << endl;

    delete a1;
    delete a2;
  }



//   //--------------------------------------------------------
//   seq1 = makeSequence( "AAAAACCCCCAAAAA");
//   seq2 = makeSequence( "AAAAAKKKKKAAAAA");
//   ali = a->Align( seq1, seq2, ali);
//   cout << *seq1 << endl << *seq2 << endl << *ali;
//   cout << "Score: " << ali->getScore() << " (should be: 37 = 10 * 4 - 2 * 1 - 10 *  0.1)" <<endl;
//   if (mode == PAIR) writePairAlignment( cout, seq1, seq2, ali);
//   delete seq1;
//   delete seq2;

//   //--------------------------------------------------------
//   seq1 = makeSequence( "AAAACAAAAAA" );
//   seq2 = makeSequence( "AAAAAAAAAA");
//   ali = a->Align( seq1, seq2, ali);
//   cout << *seq1 << endl << *seq2 << endl << *ali;
//   cout << "Score: " << ali->getScore() << " (should be: 37 = 10 * 4 - 2 * 1 - 10 *  0.1)" <<endl;
//   if (mode == PAIR) writePairAlignment( cout, seq1, seq2, ali);
//   delete seq1;
//   delete seq2;

//   //--------------------------------------------------------
//   seq1 = makeSequence( "AAAAKKKKAAAAAA" );
//   seq2 = makeSequence( "KKAKKKGGG");
//   ali = a->Align( seq1, seq2, ali);
//   cout << *seq1 << endl << *seq2 << endl << *ali;
//   cout << "Score: " << ali->getScore() << " (should be: 19 = 1 * 4 + 3 * 5)" <<endl;
//   if (mode == PAIR) writePairAlignment( cout, seq1, seq2, ali);
//   delete seq1;
//   delete seq2;

//   //--------------------------------------------------------
//   // check wraparound-alignment
//   seq1 = makeSequence( "YKKAAAAAAYKKAAAAAAYKKAAAAAAYKKAAAAAA" );
//   seq2 = makeSequence( "KKAAAAAA");
//   ali = a->Align( seq1, seq2, ali);
//   cout << *seq1 << endl << *seq2 << endl << *ali;
//   ali = a->Align( seq2, seq1, ali);
//   cout << *ali;
//   delete seq1;
//   delete seq2;

  //--------------------------------------------------------
  cout << "Test finished..." << endl;

}
 */
int main () {

	HSubstitutionMatrix matrix = makeSubstitutionMatrix( 
			getDefaultEncoder()->getAlphabetSize(), 
			10, -1);
	
	setDefaultSubstitutionMatrix( matrix );

	HAlignandum seq1 = makeSequence( "AAAAACCCCCAAAAA" );
	HAlignandum seq2 = makeSequence( "CCCCC" );
	HAlignandum seq3 = makeSequence( "CCCKCCC" );
	HAlignandum seq4 = makeSequence( "AAAAACCACCAAAAA" );
	HAlignandum seq5 = makeSequence( "KKKACACACKKK");
	HAlignandum seq6 = makeSequence( "AC");  
	HAlignandum seq7 = makeSequence( "AAAAAAACCCCAAAAAAA" );
	HAlignandum seq8 = makeProfile( "AAAAAAACCCCAAAAAAA", 1);

	Score gop = -12;
	Score gep = -2;
	
	{
		std::cout << "--- testing AlignatorDPFull (global mode) " << std::endl;
		HAlignator a = makeAlignatorDPFull( ALIGNMENT_GLOBAL, gop, gep, true, true );
		testPairwiseAlignment( 1, a, seq1, seq1, 0, 15, "+15",    0, 15,     "+15", 150 );
		testPairwiseAlignment( 2, a, seq1, seq2, 5, 10, "+5",     0,  5,      "+5", 6 );
		testPairwiseAlignment( 3, a, seq2, seq1, 0,  5, "+5",     5, 10,      "+5", 6 );
		testPairwiseAlignment( 4, a, seq2, seq3, 0,  5, "+3-2+2", 0,  7,      "+7", 34 );
		testPairwiseAlignment( 5, a, seq1, seq4, 0,  15, "+15",     0,  15,  "+15", 139 );

		HIterator2D i = makeIterator2DBanded( 0, 0);
		a->setIterator2D( i );

		testPairwiseAlignment( 6, a, seq1, seq1, 0, 15, "+15",    0, 15,     "+15", 150 );
		testPairwiseAlignment( 7, a, seq1, seq2, 0,  5, "+5",     0,  5,      "+5", -5 );
		testPairwiseAlignment( 8, a, seq2, seq3, 0,  5, "+5",     0,  5,      "+5", 39 );
		testPairwiseAlignment( 9, a, seq1, seq4, 0,  15, "+15",     0,  15,  "+15", 139 );

	}

	{
		std::cout << "--- testing AlignatorDPFull (local mode)" << std::endl;
		HAlignator a = makeAlignatorDPFull( ALIGNMENT_LOCAL, gop, gep );
		testPairwiseAlignment(11, a, seq1, seq1, 0, 15, "+15",    0, 15,     "+15", 150 );
		testPairwiseAlignment(12, a, seq1, seq2, 5, 10, "+5",     0,  5,      "+5", 50 );
		testPairwiseAlignment(13, a, seq2, seq1, 0,  5, "+5",     5, 10,      "+5", 50 );
		testPairwiseAlignment(14, a, seq2, seq3, 0,  5, "+5",     0,  5,      "+5", 39 );

		HIterator2D i = makeIterator2DBanded( 0, 0);
		a->setIterator2D( i );

		testPairwiseAlignment( 15, a, seq1, seq1, 0, 15, "+15",    0, 15,     "+15", 150 );
		testPairwiseAlignment( 16, a, seq1, seq2, NO_POS, NO_POS, "",     NO_POS,  NO_POS,      "", 0 );
		testPairwiseAlignment( 17, a, seq2, seq3, 0,  5, "+5",     0,  5,      "+5", 39 );
		testPairwiseAlignment( 18, a, seq1, seq4, 0,  15, "+15",     0,  15,  "+15", 139 );
	}
	{
		std::cout << "--- testing AlignatorDPWrap (local mode)" << std::endl;
		HAlignator a = makeAlignatorDPFull( ALIGNMENT_WRAP, gop, gep );
		testWrappedAlignment(21, a, seq5, seq6, "-7:-0+2;-5:-0+2;-3:-0+2", 60 );
		testWrappedAlignment(22, a, seq6, seq5, "3:-0+2", 20 );    
	}

	{ 
		std::cout << "--- testing setting of range ---" << std::endl;
		HAlignator a(makeAlignatorDPFull( ALIGNMENT_GLOBAL, gop, gep, true, true ));
		seq1->useSegment(5,8);
		testPairwiseAlignment( 23, a, seq1, seq2, 5, 8, "+3", 0, 3, "+3", 14 );
		seq1->useSegment(); 
	}

	{ 
		std::cout << "--- testing sequence/profile alignment ---" << std::endl;
		HAlignator a(makeAlignatorDPFull( ALIGNMENT_GLOBAL, gop, gep, true, true ));
		testPairwiseAlignment( 24, a, seq7, seq8, 0, 18, "+18", 0, 18, "+18", 18 );
	}

	{
		// here, gap costs are -10, -2
		std::cout << "--- testing AlignatorGroupies " << std::endl;
		HAlignator a = makeAlignatorGroupies();
		testPairwiseAlignment( 1, a, seq1, seq1, 0, 15, "+15",    0, 15,     "+15", 150 );
		testPairwiseAlignment( 2, a, seq1, seq2, 5, 10, "+5",     0,  5,      "+5", 50 );
		testPairwiseAlignment( 3, a, seq2, seq1, 0,  5, "+5",     5, 10,      "+5", 50 );
		testPairwiseAlignment( 4, a, seq2, seq3, 0,  5, "+2-2+3", 0,  7,      "+7", 34 );
		testPairwiseAlignment( 5, a, seq1, seq4, 0,  15, "+8-1+7",0,  15,  "+7-1+8", 112 );
	}
	
	
	{ 
		
		HSubstitutionMatrix new_matrix(makeSubstitutionMatrix( 
				getDefaultEncoder()->getAlphabetSize(), 10, -10));
		setDefaultSubstitutionMatrix( new_matrix );

		std::cout << "--- testing iterative alignment ---" << std::endl;
		HAlignator a(makeAlignatorDPFull( ALIGNMENT_LOCAL, -10.0, -2.0 ));

		HAlignment result = makeAlignmentVector();
		HAlignator alignator = makeAlignatorIterative( a, 1.0 );
		{
			HAlignandum row = makeSequence( "AAACCCCCCCCAAACCCCCCCAAACCCCCCCAAACCCCCCCAAA" );
			HAlignandum col = makeSequence( "AAAKKKKKKKKAAAKKKKKKKAAAKKKKKKKAAAKKKKKKKAAA" );
		
			alignator->align( result, row, col );		
		
			std::cout << *result << std::endl;
			std::cout << AlignmentFormatExplicit( result, row, col ) << std::endl;
		}
		{
			HAlignandum row = makeProfile( "AAACCCCCCCCAAACCCCCCCAAACCCCCCCAAACCCCCCCAAAAAACCCCCCCCAAACCCCCCCAAACCCCCCCAAACCCCCCCAAA", 2,
					getDefaultEncoder(),
					makeWeightor(),
					makeRegularizor(),
					makeLogOddor());
			HAlignandum col = makeProfile( "AAAKKKKKKKKAAAKKKKKKKAAAKKKKKKKAAAKKKKKKKAAAAAAKKKKKKKKAAAKKKKKKKAAAKKKKKKKAAAKKKKKKKAAA", 2,
					getDefaultEncoder(),
					makeWeightor(),
					makeRegularizor(),
					makeLogOddor());
		
			alignator->align( result, row, col );		
		
			std::cout << *result << std::endl;
			std::cout << AlignmentFormatExplicit( result, row, col ) << std::endl;			
		}
	}
	
	//   cout << "---------------------testing AlignatorFullDPWrap----------------------------------" << endl;
	//   a1 = makeAlignatorFullDPWrap( -1, -0.1 );
	//   ali= makeAlignmentSet();
	//   testAlignator( a1, ali, PAIR   );
	//   delete a1;
	//   delete ali;

	//    cout << "---------------------testing AlignatorFullDPGlobal: end-gap-penalties on both sides----------------------------------" << endl;
	//    a1 = makeAlignatorFullDPGlobal( -1, -0.1, true, true );
	//    ali= makeAlignmentSet();
	//    testAlignator( a1, ali, PAIR );
	//    delete a1;
	//    delete ali;

	//    cout << "---------------------testing AlignatorFullDPGlobal: no right end-gap-penalty----------------------------------" << endl;
	//    a1 = makeAlignatorFullDPGlobal( -1, -0.1, true, false );
	//    ali= makeAlignmentSet();
	//    testAlignator( a1, ali, PAIR );
	//    delete a1;
	//    delete ali;

	//    cout << "---------------------testing AlignatorFullDPGlobal: no left end gap-penalty----------------------------------" << endl;
	//    a1 = makeAlignatorFullDPGlobal( -1, -0.1, false, true );
	//    ali= makeAlignmentSet();
	//    testAlignator( a1, ali, PAIR );
	//    delete a1;
	//    delete ali;

	//    cout << "---------------------testing AlignatorFullDPGlobal: no end gap-penalty----------------------------------" << endl;
	//    a1 = makeAlignatorFullDPGlobal( -1, -0.1, false, false );
	//    ali= makeAlignmentSet();
	//    testAlignator( a1, ali, PAIR );
	//    delete a1;
	//    delete ali;


	//   cout << "---------------------testing AlignatorIdentity----------------------------------" << endl;
	//   a1 = makeAlignatorIdentity();
	//   ali= makeAlignmentMatrixRow();
	//   testAlignator( a1, ali, PAIR );
	//   delete a1;
	//   delete ali;

	//   cout << "---------------------testing AlignatorSimilarity----------------------------------" << endl;
	//   a1 = makeAlignatorSimilarity();		
	//   ali= makeAlignmentMatrixRow();
	//   testAlignator( a1, ali, MATRIX);
	//   delete a1;
	//   delete ali;

	//   cout << "---------------------testing AlignatorTuples----------------------------------" << endl;
	//   a1 = makeAlignatorTuples();		// this uses default values for ktuple-size and substitution matrix
	//   ali= makeAlignmentMatrixRow();
	//   testAlignator( a1, ali, MATRIX );
	//   delete a1;
	//   delete ali;

	/*
  {
    cout << "---------------------testing AlignatorDots----------------------------------" << endl;

    // this uses default values for ktuple-size and substitution matrix
    Alignator * t  = makeAlignatorTuples( 1 );
    Alignator * a  = makeAlignatorDotsSquared( -1, -0.1, t);		
    Alignment  * ali= makeAlignmentVector();

    testPairAlignator( a, ali );
    delete a;
    delete t;
    delete ali;
  }
	 */

}
