/*
  alignlib - a library for aligning protein sequences

  $Id: test_AlignatorDots.cpp,v 1.3 2004/06/02 12:11:38 aheger Exp $

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


#include "HelpersAlignandum.h"
#include "HelpersScorer.h"

#include "Alignator.h"
#include "HelpersAlignator.h"

#include "Alignandum.h"

#include "Alignment.h"
#include "HelpersAlignment.h"

using namespace std;

using namespace alignlib;

void test() 
{
  {
        cout << "---------------------Testing Alignator----------------------------------" << endl;
        // build alignment
    
    HAlignment dotplot( makeAlignmentMatrixRow() );
    
    AlignmentFormatDiagonals f;
    f.mAlignment = std::string("-321:-89+10;-320:-78+11-11+13;-315:-75+3;-311:-70+5;-310:-47+23;-309:-40+7;-305:-28+12;-304:-27+1;-303:-13+5;-302:-4+9;-301:-20+7;-169:-98+2;-167:-95+3-4+8;-166:-94+1;-164:-60+15;-163:-76+18;-156:-35+25;-155:-5+10;-154:-16+19;-24:-53+2;-23:-92+7;-22:-31+4-57+5;-21:-31+81;-20:-30+6-23+1;-19:-30+6;-18:-33+3;-17:-0+1-32+3;-16:-0+36;-15:-37+1");
    f.mRowFrom = 0;
    f.mColFrom = 0;
    f.copy( dotplot );
    
	HAlignandum seq2 (makeSequence( "VVRCDIQMTQSPASLSASVGETVTITCGASENIYGALNWYQRKQGKSPQLLIYGATNLADGMSSRFSGSGSGRQYSLKISSLHPDDVATYYCQNVLSAPWTFGGGTKLEIKRAD",
			getDefaultEncoder()));
	HAlignandum seq1 (makeSequence( "MSVPTQVLGLLLLWLTDARCDIQMTQSPSTLSASVGDRVTITCRSSKSLLHSNGDTFLYWFQQKPGKAPKLLMYRMSNLASGVPSRFSGSGSGTEFTLTISSLQPDDFATYYCMQHLEYPFTFGQGTKVEVKRTGGGGSGGGGSGGGGSGGGGSGGGGSQIQLVQSGAEVKKPGSSVKVSCKASGYTFTDYYINWMRQAPGQGLEWIGWIDPGSGNTKYNEKFKGRATLTVDTSTNTAYMELSSLRSEDTAFYFCAREKTTYYYAMDYWGQGTLVTVSSASTKGPTSDKTHTCPPCPAPELLGGPSVFLFPPKPKDTLMISRTPEVTCVVVDVSHEDPEVKFNWYVDGVEVHNAKTKPREEQYNSTYRVVSVLTVLHQDWLNGKEYKCKVSNKALPAPIEKTISKAKGQPREPQVYTLPPSREEMTKNQVSLTCLVKGFYPSDIAVEWESNGQPENNYKTTPPVLDSDGSFFLYSKLTVDKSRWQQGNVFSCSVMHEALHNHYTQKSLSLSPGKLDPKLCYLLDGILFIYGVILTALFLRVKFSRSADAPAYQQGQNQLYNELNLGRREEYDVLDKRRGRDPEMGGKPRRKNPQEGLYNELQKDKMAEAYSEIGMKGERRRGKGHDGLYQGLSTATKDTYDALHMQALPPRRSKRSRLLHSDYMNMTPRRPGPTRKHYQPYAPPRDFAAYRS",
			getDefaultEncoder()));

  	rescoreAlignment( dotplot, seq1, seq2, makeScorer( seq1, seq2 ) );
  	// std::cout << *dotplot << endl;
	
	// make a Dottor
	HAlignator dottor(makeAlignatorPrebuilt( dotplot ));
	// make the alignator
	
	HAlignator alignator(makeAlignatorDotsSquared( dottor, -10, -2 ));
	
	// setup up map_source2dest; m2 will be in col
	HAlignment map_source2dest(makeAlignmentVector());

	alignator->align( map_source2dest, seq1, seq2 );
	
	// std::cout << *map_source2dest << endl;
	
	std::cout << AlignmentFormatExplicit( map_source2dest, seq1, seq2 ) << std::endl;
  }
  
}



int main () {

  test();
//    {
//          cout << "---------------------Testing Alignator----------------------------------" << endl;
//          // build alignment
//          alignlib::Alignment * dotplot = fillAlignmentIdentity( makeAlignmentMatrixRow(), 10, 20, 0);
//  	rescoreAlignment( dotplot, 1);
//  	cout << *dotplot << endl;
  
//  	Alignandum * pro1 = makeProfile( "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA", 1);
//  	Alignandum * pro2 = makeProfile( "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA", 1);
	
//  	pro1->useSegment(10,20);
//  	pro2->useSegment(10,20);
	
//  	// make a Dottor
//  	alignlib::Alignator * dottor = makeAlignatorPublishAlignment( dotplot );
	
//  	// make the alignator
//  	alignlib::Alignator * alignator = makeAlignatorDotsSquared( -4, -0.4, dottor );
	
//  	// setup up map_source2dest; m2 will be in col
//  	Alignment *map_source2dest = makeAlignmentSet();
	
//  	alignator->Align( pro1, pro2, map_source2dest);
	
//  	cout << *map_source2dest << endl;
	
//  	delete alignator;
//  	delete pro1;
//  	delete pro2;
//  	delete dottor;
//  	delete dotplot;
//     }


    
//     {
// 	cout << "-------------------------------Test 2-----------------------------------------" << endl;
// 	Alignandum * pro1 = makeProfile( "AAGGGGGGGGGGG", 1);
// 	Alignandum * pro2 = makeProfile( "AAGGGGGGGGGGGYYAAGGGGGGGGGGGYYAAGGGGGGGGGGGYYAAGGGGGGGGGGG", 1);
// 	Alignator * dottor = makeAlignatorIdentity();
// 	Alignator * alignator = makeAlignatorDotsWrap( -4, -0.1, dottor);
// 	Alignment * ali = makeAlignmentSet();
// 	cout << "pro1 -> pro2" << endl;
// 	alignator->Align( pro1, pro2, ali);
// 	cout << *ali << endl;
// 	cout << "pro2 -> pro1" << endl;
// 	alignator->Align( pro2, pro1, ali);
// 	cout << *ali << endl;

// 	delete alignator;
// 	delete pro1;
// 	delete pro2;
// 	delete dottor;
	
//     }

//      {
//  	cout << "-------------------------------Test 2-----------------------------------------" << endl;
//  // 	Alignandum * seq1 = makeSequence( "XXGXXGX" );
//  // 	Alignandum * seq2 = makeSequence( "XXGXXGX" );
//  	Alignment * a = makeAlignmentMatrixRow();

//   	Alignandum * seq1 = makeSequence( "XXGXXGXXXXAXXAX" );
//   	Alignandum * seq2 = makeSequence( "XXGXXGXAXXAXXX" );
//  	cout << "seq1: " << *seq1 << endl;
//  	cout << "seq2: " << *seq2 << endl;
//  	Alignator * dottor = makeAlignatorIdentity();
//  	dottor->Align( seq1, seq2, a);
//  	cout << *a << endl;

//  	Alignator * alignator = makeAlignatorDotsDiagonals( 0, -0.1, dottor);
//  	Alignment * ali = makeAlignmentSet();
//  	cout << "seq1 -> seq2" << endl;
//  	alignator->Align( seq1, seq2, ali);
//  	cout << *ali << endl;
//  	cout << "seq2 -> seq1" << endl;
//  	alignator->Align( seq2, seq1, ali);
//  	cout << *ali << endl;

//  	delete alignator;
//  	delete seq1;
//  	delete seq2;
//  	delete dottor;

//      }

  return EXIT_SUCCESS;
}
