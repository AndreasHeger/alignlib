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

#define BOOST_TEST_MODULE
#include <boost/test/included/unit_test.hpp>
using boost::unit_test::test_suite;

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

	// std::cout << "======== test " << test_id << " =========" << std::endl;

	// std::cout << *benchmark_row << std::endl;
	// std::cout << *benchmark_col << std::endl;

	HAlignment result(makeAlignmentVector());

	a->align( result, benchmark_row, benchmark_col );

	std::string r_row(row);
	std::string r_col(col);

	// std::cout << AlignmentFormatExplicit( result, benchmark_row, benchmark_col ) << std::endl ;

	// std::cout << "result=" << *result << std::endl;

	AlignmentFormatEmissions format( result );

	BOOST_CHECK_EQUAL(result->getScore(), score);
	BOOST_CHECK_EQUAL(row_from, result->getRowFrom());
	BOOST_CHECK_EQUAL(row_to, result->getRowTo());
	BOOST_CHECK_EQUAL(col_from, result->getColFrom());
	BOOST_CHECK_EQUAL(col_to, result->getColTo());
	BOOST_CHECK_EQUAL(format.mRowAlignment, r_row);
	BOOST_CHECK_EQUAL(format.mColAlignment, r_col );

	if ( result->getScore() == score &&
			row_from == result->getRowFrom() &&
			row_to == result->getRowTo() &&
			col_from == result->getColFrom() &&
			col_to == result->getColTo() &&
			format.mRowAlignment == r_row &&
			format.mColAlignment == r_col )
	{
		// std::cout << "test " << test_id << " success" << std::endl;
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

	// std::cout << "======== test " << test_id << " =========" << std::endl;
	// std::cout << *benchmark_row << std::endl;
	// std::cout << *benchmark_col << std::endl;

	HAlignment result(makeAlignmentMatrixDiagonal());

	a->align( result, benchmark_row, benchmark_col );

	std::string r_row(r_ali);

	AlignmentFormatDiagonals format( result );

	if ( result->getScore() == score &&
			format.mAlignment == r_ali )
	{
		// std::cout << "test " << test_id << " success" << std::endl;
		return true;
	}
	else
	{
		std::cout << *result << std::endl;
		std::cout << r_ali
		<< (( format.mAlignment == r_ali ) ? " == " : " != ")
		<< format.mAlignment
		<< std::endl;
		std::cout << "test " << test_id << " failure" << std::endl;

		return false;
	}
}

BOOST_AUTO_TEST_CASE( global_alignment2)
{
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
		// std::cout << "--- testing AlignatorDPFull (global mode) " << std::endl;
		HAlignator a = makeAlignatorDPFull( ALIGNMENT_GLOBAL, gop, gep, true, true, true, true );
		testPairwiseAlignment( 1, a, seq1, seq1, 0, 15, "+15",    0, 15,     "+15", 150 );
		testPairwiseAlignment( 2, a, seq1, seq2, 5, 10, "+5",     0,  5,      "+5", 6 );
		testPairwiseAlignment( 3, a, seq2, seq1, 0,  5, "+5",     5, 10,      "+5", 6 );
		testPairwiseAlignment( 4, a, seq2, seq3, 0,  5, "+3-2+2", 0,  7,      "+7", 34 );
		testPairwiseAlignment( 5, a, seq1, seq4, 0,  15, "+15",     0,  15,  "+15", 139 );

		HIterator2D i = makeIterator2DBanded( 0, 0);
		a->cloneToolkit();
		a->getToolkit()->setIterator2D( i );

		testPairwiseAlignment( 6, a, seq1, seq1, 0, 15, "+15",    0, 15,     "+15", 150 );
		testPairwiseAlignment( 7, a, seq1, seq2, 0,  5, "+5",     0,  5,      "+5", -5 );
		testPairwiseAlignment( 8, a, seq2, seq3, 0,  5, "+5",     0,  5,      "+5", 39 );
		testPairwiseAlignment( 9, a, seq1, seq4, 0,  15, "+15",     0,  15,  "+15", 139 );

	}


	{
		// std::cout << "--- testing setting of range ---" << std::endl;
		HAlignator a(makeAlignatorDPFull( ALIGNMENT_GLOBAL, gop, gep, true, true ));
		seq1->useSegment(5,8);
		testPairwiseAlignment( 23, a, seq1, seq2, 5, 8, "+3", 0, 3, "+3", 14 );
		seq1->useSegment();
	}

	{
		// std::cout << "--- testing sequence/profile alignment ---" << std::endl;
		HAlignator a(makeAlignatorDPFull( ALIGNMENT_GLOBAL, gop, gep, true, true ));
		testPairwiseAlignment( 24, a, seq7, seq8, 0, 18, "+18", 0, 18, "+18", 18 );
	}

}


BOOST_AUTO_TEST_CASE( alignment_wrap)
{
	HSubstitutionMatrix matrix = makeSubstitutionMatrix(
			getDefaultEncoder()->getAlphabetSize(),
			10, -1);

	setDefaultSubstitutionMatrix( matrix );

	HAlignandum seq1 = makeSequence( "KKKACACACKKK");
	HAlignandum seq2 = makeSequence( "AC");

	Score gop = -12;
	Score gep = -2;

	// std::cout << "--- testing AlignatorDPWrap (local mode)" << std::endl;
	HAlignator a = makeAlignatorDPFull( ALIGNMENT_WRAP, gop, gep );
	testWrappedAlignment(21, a, seq1, seq2, "-7:-0+2;-5:-0+2;-3:-0+2", 60 );
	testWrappedAlignment(22, a, seq2, seq1, "3:-0+2", 20 );
}

BOOST_AUTO_TEST_CASE( alignment_wrap2)
{
	HSubstitutionMatrix matrix = makeSubstitutionMatrix(
			getDefaultEncoder()->getAlphabetSize(),
			10, -1);

	setDefaultSubstitutionMatrix( matrix );

	HAlignandum seq1 = makeSequence( "KKKACACACKKK");
	HAlignandum seq2 = makeSequence( "PACP");

	Score gop = -12;
	Score gep = -2;

	// std::cout << "--- testing AlignatorDPWrap (local mode)" << std::endl;
	HAlignator a = makeAlignatorDPFull( ALIGNMENT_WRAP, gop, gep );
	testWrappedAlignment(23, a, seq1, seq2, "-6:-0+3;-2:-1+3", 38 );
	testWrappedAlignment(24, a, seq2, seq1, "2:-1+2", 20 );
}

BOOST_AUTO_TEST_CASE( iterative_alignment )
{
	HSubstitutionMatrix new_matrix(makeSubstitutionMatrix(
			getDefaultEncoder()->getAlphabetSize(), 10, -10));
	setDefaultSubstitutionMatrix( new_matrix );

	// std::cout << "--- testing iterative alignment ---" << std::endl;
	HAlignator a(makeAlignatorDPFull( ALIGNMENT_LOCAL, -10.0, -2.0 ));

	HAlignment result = makeAlignmentVector();
	HAlignator alignator = makeAlignatorIterative( a, 1.0 );
	{
		HAlignandum row = makeSequence( "AAACCCCCCCCAAACCCCCCCAAACCCCCCCAAACCCCCCCAAA" );
		HAlignandum col = makeSequence( "AAAKKKKKKKKAAAKKKKKKKAAAKKKKKKKAAAKKKKKKKAAA" );

		alignator->align( result, row, col );

		// std::cout << *result << std::endl;
		// std::cout << AlignmentFormatExplicit( result, row, col ) << std::endl;
	}
	{
		HAlignandum row = makeProfile( "AAACCCCCCCCAAACCCCCCCAAACCCCCCCAAACCCCCCCAAAAAACCCCCCCCAAACCCCCCCAAACCCCCCCAAACCCCCCCAAA", 2 );
		HAlignandum col = makeProfile( "AAAKKKKKKKKAAAKKKKKKKAAAKKKKKKKAAAKKKKKKKAAAAAAKKKKKKKKAAAKKKKKKKAAAKKKKKKKAAAKKKKKKKAAA", 2 );

		alignator->align( result, row, col );

		// std::cout << *result << std::endl;
		// std::cout << AlignmentFormatExplicit( result, row, col ) << std::endl;
	}
}

BOOST_AUTO_TEST_CASE( groupies_alignment)
{
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
		// here, gap costs are -10, -2
		// std::cout << "--- testing AlignatorGroupies " << std::endl;
		HAlignator a = makeAlignatorGroupies();
		testPairwiseAlignment( 31, a, seq1, seq1, 0, 15, "+15",    0, 15,     "+15", 150 );
		//testPairwiseAlignment( 32, a, seq1, seq2, 5, 10, "+5",     0,  5,      "+5", 50 );
		//testPairwiseAlignment( 33, a, seq2, seq1, 0,  5, "+5",     5, 10,      "+5", 50 );
		//testPairwiseAlignment( 34, a, seq2, seq3, 0,  5, "+2-2+3", 0,  7,      "+7", 36 );
		//testPairwiseAlignment( 35, a, seq1, seq4, 0,  15, "+15",0,  15,  "+15", 139 );
	}
}

BOOST_AUTO_TEST_CASE( local_alignment)
{
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
		// std::cout << "--- testing AlignatorDPFull (local mode)" << std::endl;
		HAlignator a = makeAlignatorDPFull( ALIGNMENT_LOCAL, gop, gep );
		testPairwiseAlignment(11, a, seq1, seq1, 0, 15, "+15",    0, 15,     "+15", 150 );
		testPairwiseAlignment(12, a, seq1, seq2, 5, 10, "+5",     0,  5,      "+5", 50 );
		testPairwiseAlignment(13, a, seq2, seq1, 0,  5, "+5",     5, 10,      "+5", 50 );
		testPairwiseAlignment(14, a, seq2, seq3, 0,  5, "+5",     0,  5,      "+5", 39 );
		HIterator2D i = makeIterator2DBanded( 0, 0);
		a->cloneToolkit();
		a->getToolkit()->setIterator2D( i );

		testPairwiseAlignment( 15, a, seq1, seq1, 0, 15, "+15",    0, 15,     "+15", 150 );
		testPairwiseAlignment( 16, a, seq1, seq2, NO_POS, NO_POS, "",     NO_POS,  NO_POS,      "", 0 );
		testPairwiseAlignment( 17, a, seq2, seq3, 0,  5, "+5",     0,  5,      "+5", 39 );
		testPairwiseAlignment( 18, a, seq1, seq4, 0,  15, "+15",     0,  15,  "+15", 139 );
	}

}


BOOST_AUTO_TEST_CASE( global_alignment)
{

	{
		HSubstitutionMatrix matrix = makeSubstitutionMatrix(
				getDefaultEncoder()->getAlphabetSize(),
				1, -1);

		setDefaultSubstitutionMatrix( matrix );

		HAlignandum seq1 = makeSequence( "AAAAACCCCCAAAAACCCCCAAAAA" );
		HAlignandum seq2 = makeSequence( "CCCCCAAAAA" );

		Score gop = -2;
		Score gep = -1;

		{
			HAlignator a = makeAlignatorDPFull( ALIGNMENT_GLOBAL, gop, gep, false, false, false, false );
			testPairwiseAlignment( 51, a, seq1, seq2, 5,  15, "+10",     0,  10,  "+10", 10 );
			testPairwiseAlignment( 52, a, seq2, seq1, 0,  10, "+10",     5,  15,  "+10", 10);
		}
		{
			HAlignator a = makeAlignatorDPFull( ALIGNMENT_GLOBAL, gop, gep, true, false, false, false );
			testPairwiseAlignment( 53, a, seq1, seq2, 5,  15, "+10",     0,  10,  "+10", 10 );
			testPairwiseAlignment( 54, a, seq2, seq1, 5,  10, "+5",     0,  5,  "+5", 5 );
		}
		{
			HAlignator a = makeAlignatorDPFull( ALIGNMENT_GLOBAL, gop, gep, false, true, false, false );
			testPairwiseAlignment( 55, a, seq1, seq2, 5,  15, "+10",     0,  10,  "+10", 10 );
			testPairwiseAlignment( 56, a, seq2, seq1, 0,  10, "+10",    15,  25,  "+10", 10 );
		}
		{
			HAlignator a = makeAlignatorDPFull( ALIGNMENT_GLOBAL, gop, gep, true, true, false, false );
			testPairwiseAlignment( 57, a, seq1, seq2, 5,  15, "+10",     0,  10,  "+10", 10 );
			testPairwiseAlignment( 58, a, seq2, seq1, 0,  10, "+10",    15,  25,  "+10", -7 );
		}
		{
			HAlignator a = makeAlignatorDPFull( ALIGNMENT_GLOBAL, gop, gep, false, false, true, false );
			testPairwiseAlignment( 59, a, seq2, seq1, 0,  10, "+10",     5,  15,  "+10", 10 );
			testPairwiseAlignment( 60, a, seq1, seq2, 0,  5,  "+5",      5,  10,  "+5", 5 );
		}
		{
			HAlignator a = makeAlignatorDPFull( ALIGNMENT_GLOBAL, gop, gep, false, false, false, true );
			testPairwiseAlignment( 61, a, seq2, seq1, 0,  10, "+10",     5,  15,  "+10", 10 );
			testPairwiseAlignment( 62, a, seq1, seq2, 15,  25, "+10",    0,  10,  "+10", 10 );
		}
		{
			HAlignator a = makeAlignatorDPFull( ALIGNMENT_GLOBAL, gop, gep, false, false, true, true );
			testPairwiseAlignment( 63, a, seq2, seq1, 0,  10, "+10",     5,  15,  "+10", 10 );
			testPairwiseAlignment( 64, a, seq1, seq2,15,  25, "+10",     0,  10,  "+10", -7 );
		}
		{
			HAlignator a = makeAlignatorDPFull( ALIGNMENT_GLOBAL, gop, gep, true, true, true, true );
			testPairwiseAlignment( 65, a, seq2, seq1,  0,  10, "+10",   15,  25, "+10",      -7 );
			testPairwiseAlignment( 66, a, seq1, seq2, 15,  25, "+10",    0,  10, "+10",      -7 );
		}
	}

	{
		HAlignandum seq1 = makeSequence( "AAAAACCCCCAAAAA" );
		HAlignandum seq2 = makeSequence( "YYYYYCCCCCYYYYY" );

		Score gop = -2;
		Score gep = -1;

		{
			HAlignator a = makeAlignatorDPFull( ALIGNMENT_GLOBAL, gop, gep, false, false, false, false );
			testPairwiseAlignment( 70, a, seq1, seq2, 5,  10, "+5",     5,  10,  "+5", 5 );
			testPairwiseAlignment( 71, a, seq2, seq1, 5,  10, "+5",     5,  10,  "+5", 5 );
		}

	}
}

BOOST_AUTO_TEST_CASE( alignment_backtranslation )
{
	HSubstitutionMatrix matrix = makeSubstitutionMatrixBackTranslation(
			2, -10, 1, getDefaultEncoder() );

	setDefaultSubstitutionMatrix( matrix );

	HAlignandum seq1 = makeSequence( "CAYACWGCWTGYCCWWGWWTWGTWGTWWT" );
	HAlignandum seq2 = makeSequence( "CATACAGCGTGCCCTGCGGCTGGTGGTGCT" );

	Score gop = -2;
	Score gep = -1;

	{
		HAlignator a = makeAlignatorDPFull( ALIGNMENT_GLOBAL, gop, gep, false, false, false, false );
		testPairwiseAlignment( 81, a, seq1, seq2, 0,  29, "+16-1+13", 0,  30,  "+30", 45 );
		testPairwiseAlignment( 82, a, seq2, seq1, 0,  30, "+30",      0,  29,  "+16-1+13", 45);
	}
}

BOOST_AUTO_TEST_CASE( alignment_backtranslation_with_gaps )
{
        HSubstitutionMatrix matrix = makeSubstitutionMatrixBackTranslation(
                        2, -10, 1, getDefaultEncoder() );

        setDefaultSubstitutionMatrix( matrix );

        HAlignandum seq1 = makeSequence( "TTGGCTGTGAATCTATTCCATCTGTGAGTG" );
        HAlignandum seq2 = makeSequence( "WTWGCWGTWWWW---ATWCCWWWWGTWWWW" );

        Score gop = -2;
        Score gep = -1;

        {
        	HAlignator a = makeAlignatorDPFull( ALIGNMENT_GLOBAL, gop, gep, false, false, false, false );
        	testPairwiseAlignment( 83, a, seq1, seq2, 0,  29, "+14-3+15", 0,  30,  "+12-2+18", 29 );
        	testPairwiseAlignment( 84, a, seq2, seq1, 0,  30, "+15-2+15", 0,  29,  "+12-3+17", 29 );
        }
}


