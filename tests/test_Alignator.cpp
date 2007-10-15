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

#include "SubstitutionMatrix.h"
#include "HelpersSubstitutionMatrix.h"

#include "HelpersSequence.h"

#include "HelpersProfile.h"

#include "Alignator.h"
#include "HelpersAlignator.h"

#include "HelpersIterator2D.h"

#include "Alignata.h"
#include "HelpersAlignata.h"

#ifdef WITH_DMALLOC
#include <dmalloc.h>
#endif

using namespace std;

#ifdef WITH_DMALLOC
#include <dmalloc.h>
#endif

using namespace alignlib;

// typedef enum MODE { PAIR, MATRIX };



bool TestPairwiseAlignment(int test_id,
			   alignlib::Alignator * a, 
			   const Alignandum * benchmark_row,
			   const Alignandum * benchmark_col,
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
  
  alignlib::Alignata * result = makeAlignataVector();

  a->Align( benchmark_row, benchmark_col, result);
  
  std::string b_row;
  std::string b_col;
  std::string r_row(row);
  std::string r_col(col);

  writePairAlignment( std::cout, benchmark_row, benchmark_col, result );

  writeAlignataCompressed( result, b_row, b_col );
  
  if ( result->getScore() == score &&
       b_row == r_row &&
       b_col == r_col )
    {
      std::cout << "test " << test_id << " success" << std::endl;
      return true;
    }
  else
    {
      std::cout << row_from << "-" << row_to << ":" << r_row
		<< (( b_row == r_row ) ? " == " : " != ")
		<< result->getRowFrom() << "-" << result->getRowTo() << ":" << b_row << "\t"
		<< col_from << "-" << col_to << ":" << r_col
		<< (( b_col == r_col ) ? " == " : " != ") 
      		<< result->getColFrom() << "-" << result->getColTo() << ":" << b_col << "\t"
		<< score << (( score == result->getScore() ) ? " == " : " != ") << result->getScore()
		<< std::endl;
      std::cout << "test " << test_id << " failure" << std::endl;      

      return false;
    }
  delete result;
}

bool TestWrappedAlignment(int test_id,
			  alignlib::Alignator * a, 
			  const Alignandum * benchmark_row,
			  const Alignandum * benchmark_col,
			  const char * r_ali,
			  Score score )
{

  std::cout << "======== test " << test_id << " =========" << std::endl;
  std::cout << *benchmark_row << std::endl;
  std::cout << *benchmark_col << std::endl;
  
  alignlib::Alignata * result = makeAlignataMatrixDiagonal();

  a->Align( benchmark_row, benchmark_col, result);
  
  std::string b_ali;
  std::string r_row(r_ali);

  writeAlignataCompressedDiagonal( result, b_ali );
  
  if ( result->getScore() == score &&
       b_ali == r_ali )
    {
      std::cout << "test " << test_id << " success" << std::endl;
      return true;
    }
  else
    {
      std::cout << r_ali
		<< (( b_ali == r_ali ) ? " == " : " != ")
		<< b_ali
		<< std::endl;
      std::cout << "test " << test_id << " failure" << std::endl;      

      return false;
    }
  delete result;
}


/* test pairwise alignments. For the calculation of the score it is assumed, 
   that the identity matrix is used (match = 1, mismatch = -1, gop = -0.2, gep = -0.1).
*/
/*
void TestPairAlignator( Alignator * a,
			Alignata * ali)
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

  SubstitutionMatrix * matrix = makeSubstitutionMatrixAAIdentity( 10, -1);
  
  delete setDefaultSubstitutionMatrix( matrix );

  Alignandum * seq1 = makeSequence( "AAAAACCCCCAAAAA");
  Alignandum * seq2 = makeSequence( "CCCCC");
  Alignandum * seq3 = makeSequence( "CCCKCCC");
  Alignandum * seq4 = makeSequence( "AAAAACCACCAAAAA");
  Alignandum * seq5 = makeSequence( "KKKACACACKKK");
  Alignandum * seq6 = makeSequence( "AC");  

  Score gop = -12;
  Score gep = -2;
  
  {
    std::cout << "--- testing AlignatorDPFull (global mode) " << std::endl;
    Alignator * a = makeAlignatorDPFull( ALIGNMENT_GLOBAL, gop, gep, true, true );
    TestPairwiseAlignment( 1, a, seq1, seq1, 1, 15, "+15",    1, 15,     "+15", 150 );
    TestPairwiseAlignment( 2, a, seq1, seq2, 6, 10, "+5",     1,  5,      "+5", 6 );
    TestPairwiseAlignment( 3, a, seq1, seq2, 1,  5, "+5",     6, 10,      "+5", 6 );
    TestPairwiseAlignment( 4, a, seq2, seq3, 1,  5, "+3-2+2", 1,  7,      "+7", 34 );
    TestPairwiseAlignment( 5, a, seq1, seq4, 1,  15, "+15",     1,  15,  "+15", 139 );

    Iterator2D * i = makeIterator2DBanded( NULL, NULL, 0, 0);
    a->setIterator2D( i );

    TestPairwiseAlignment( 6, a, seq1, seq1, 1, 15, "+15",    1, 15,     "+15", 150 );
    TestPairwiseAlignment( 7, a, seq1, seq2, 1,  5, "+5",     1,  5,      "+5", -5 );
    TestPairwiseAlignment( 8, a, seq2, seq3, 1,  5, "+5",     1,  5,      "+5", 39 );
    TestPairwiseAlignment( 9, a, seq1, seq4, 1,  15, "+15",     1,  15,  "+15", 139 );
    
    delete i;
    delete a;
  }

  {
    std::cout << "--- testing AlignatorDPFull (local mode)" << std::endl;
    Alignator * a = makeAlignatorDPFull( ALIGNMENT_LOCAL, gop, gep );
    TestPairwiseAlignment(11, a, seq1, seq1, 1, 15, "+15",    1, 15,     "+15", 150 );
    TestPairwiseAlignment(12, a, seq1, seq2, 6, 10, "+5",     1,  5,      "+5", 50 );
    TestPairwiseAlignment(13, a, seq1, seq2, 1,  5, "+5",     6, 10,      "+5", 50 );
    TestPairwiseAlignment(14, a, seq2, seq3, 1,  5, "+5",     1,  5,      "+5", 39 );

    Iterator2D * i = makeIterator2DBanded( NULL, NULL, 0, 0);
    a->setIterator2D( i );

    TestPairwiseAlignment( 15, a, seq1, seq1, 1, 15, "+15",    1, 15,     "+15", 150 );
    TestPairwiseAlignment( 16, a, seq1, seq2, 0,  0, "",     0,  0,      "", 0 );
    TestPairwiseAlignment( 17, a, seq2, seq3, 1,  5, "+5",     1,  5,      "+5", 39 );
    TestPairwiseAlignment( 18, a, seq1, seq4, 1,  15, "+15",     1,  15,  "+15", 139 );
    
    delete i;
    delete a;
  }

  {
    std::cout << "--- testing AlignatorDPWrap (local mode)" << std::endl;
    Alignator * a = makeAlignatorDPFull( ALIGNMENT_WRAP, gop, gep );
    TestWrappedAlignment(21, a, seq5, seq6, "-7:-0+2;-5:-0+2;-3:-0+2", 60 );
    TestWrappedAlignment(22, a, seq6, seq5, "3:-0+2", 20 );    
    delete a;
  }



  
//   cout << "---------------------Testing AlignatorFullDPWrap----------------------------------" << endl;
//   a1 = makeAlignatorFullDPWrap( -1, -0.1 );
//   ali= makeAlignataSet();
//   TestAlignator( a1, ali, PAIR   );
//   delete a1;
//   delete ali;

//    cout << "---------------------Testing AlignatorFullDPGlobal: end-gap-penalties on both sides----------------------------------" << endl;
//    a1 = makeAlignatorFullDPGlobal( -1, -0.1, true, true );
//    ali= makeAlignataSet();
//    TestAlignator( a1, ali, PAIR );
//    delete a1;
//    delete ali;

//    cout << "---------------------Testing AlignatorFullDPGlobal: no right end-gap-penalty----------------------------------" << endl;
//    a1 = makeAlignatorFullDPGlobal( -1, -0.1, true, false );
//    ali= makeAlignataSet();
//    TestAlignator( a1, ali, PAIR );
//    delete a1;
//    delete ali;

//    cout << "---------------------Testing AlignatorFullDPGlobal: no left end gap-penalty----------------------------------" << endl;
//    a1 = makeAlignatorFullDPGlobal( -1, -0.1, false, true );
//    ali= makeAlignataSet();
//    TestAlignator( a1, ali, PAIR );
//    delete a1;
//    delete ali;

//    cout << "---------------------Testing AlignatorFullDPGlobal: no end gap-penalty----------------------------------" << endl;
//    a1 = makeAlignatorFullDPGlobal( -1, -0.1, false, false );
//    ali= makeAlignataSet();
//    TestAlignator( a1, ali, PAIR );
//    delete a1;
//    delete ali;


//   cout << "---------------------Testing AlignatorIdentity----------------------------------" << endl;
//   a1 = makeAlignatorIdentity();
//   ali= makeAlignataMatrixRow();
//   TestAlignator( a1, ali, PAIR );
//   delete a1;
//   delete ali;
  
//   cout << "---------------------Testing AlignatorSimilarity----------------------------------" << endl;
//   a1 = makeAlignatorSimilarity();		
//   ali= makeAlignataMatrixRow();
//   TestAlignator( a1, ali, MATRIX);
//   delete a1;
//   delete ali;

//   cout << "---------------------Testing AlignatorTuples----------------------------------" << endl;
//   a1 = makeAlignatorTuples();		// this uses default values for ktuple-size and substitution matrix
//   ali= makeAlignataMatrixRow();
//   TestAlignator( a1, ali, MATRIX );
//   delete a1;
//   delete ali;

  /*
  {
    cout << "---------------------Testing AlignatorDots----------------------------------" << endl;

    // this uses default values for ktuple-size and substitution matrix
    Alignator * t  = makeAlignatorTuples( 1 );
    Alignator * a  = makeAlignatorDotsSquared( -1, -0.1, t);		
    Alignata  * ali= makeAlignataVector();

    TestPairAlignator( a, ali );
    delete a;
    delete t;
    delete ali;
  }
  */

  delete seq1;
  delete seq2;
  delete seq3;
  delete seq4;
  
  delete matrix;
  
}
