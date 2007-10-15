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


#include "HelpersSequence.h"

#include "HelpersProfile.h"

#include "Alignator.h"
#include "HelpersAlignator.h"

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

void test() {
  {
        cout << "---------------------Testing Alignator----------------------------------" << endl;
        // build alignment

        alignlib::Alignata * dotplot = fillAlignataCompressedDiagonal( makeAlignataMatrixRow(), 
								       "-321:-89+10;-320:-78+11-11+13;-315:-75+3;-311:-70+5;-310:-47+23;-309:-40+7;-305:-28+12;-304:-27+1;-303:-13+5;-302:-4+9;-301:-20+7;-169:-98+2;-167:-95+3-4+8;-166:-94+1;-164:-60+15;-163:-76+18;-156:-35+25;-155:-5+10;-154:-16+19;-24:-53+2;-23:-92+7;-22:-31+4-57+5;-21:-31+81;-20:-30+6-23+1;-19:-30+6;-18:-33+3;-17:-0+1-32+3;-16:-0+36;-15:-37+1"
								       );
  
	Alignandum * seq1 = makeSequence( "MSVPTQVLGLLLLWLTDARCDIQMTQSPSTLSASVGDRVTITCRSSKSLLHSNGDTFLYWFQQKPGKAPKLLMYRMSNLASGVPSRFSGSGSGTEFTLTISSLQPDDFATYYCMQHLEYPFTFGQGTKVEVKRTGGGGSGGGGSGGGGSGGGGSGGGGSQIQLVQSGAEVKKPGSSVKVSCKASGYTFTDYYINWMRQAPGQGLEWIGWIDPGSGNTKYNEKFKGRATLTVDTSTNTAYMELSSLRSEDTAFYFCAREKTTYYYAMDYWGQGTLVTVSSASTKGPTSDKTHTCPPCPAPELLGGPSVFLFPPKPKDTLMISRTPEVTCVVVDVSHEDPEVKFNWYVDGVEVHNAKTKPREEQYNSTYRVVSVLTVLHQDWLNGKEYKCKVSNKALPAPIEKTISKAKGQPREPQVYTLPPSREEMTKNQVSLTCLVKGFYPSDIAVEWESNGQPENNYKTTPPVLDSDGSFFLYSKLTVDKSRWQQGNVFSCSVMHEALHNHYTQKSLSLSPGKLDPKLCYLLDGILFIYGVILTALFLRVKFSRSADAPAYQQGQNQLYNELNLGRREEYDVLDKRRGRDPEMGGKPRRKNPQEGLYNELQKDKMAEAYSEIGMKGERRRGKGHDGLYQGLSTATKDTYDALHMQALPPRRSKRSRLLHSDYMNMTPRRPGPTRKHYQPYAPPRDFAAYRS");
	Alignandum * seq2 = makeSequence( "VVRCDIQMTQSPASLSASVGETVTITCGASENIYGALNWYQRKQGKSPQLLIYGATNLADGMSSRFSGSGSGRQYSLKISSLHPDDVATYYCQNVLSAPWTFGGGTKLEIKRAD");

  	rescoreAlignment( dotplot, seq1, seq2 );
  	// cout << *dotplot << endl;
	
	// make a Dottor
	alignlib::Alignator * dottor = makeAlignatorDummy( dotplot );
	
	// make the alignator
	alignlib::Alignator * alignator = makeAlignatorDotsSquared( -10, -2, dottor );
	
	// setup up map_source2dest; m2 will be in col
	Alignata *map_source2dest = makeAlignataVector();

	for (int i = 0; i < 1000; i++) {
	  alignator->Align( seq1, seq2, map_source2dest);
	}
	
	// cout << *map_source2dest << endl;
	
	writePairAlignment( cout, seq1, seq2, map_source2dest );

	delete alignator;
	delete seq1;
	delete seq2;
	delete dottor;
	delete dotplot;
	delete map_source2dest;
  }
  
}



int main () {

  test();
//    {
//          cout << "---------------------Testing Alignator----------------------------------" << endl;
//          // build alignment
//          alignlib::Alignata * dotplot = fillAlignataIdentity( makeAlignataMatrixRow(), 10, 20, 0);
//  	rescoreAlignment( dotplot, 1);
//  	cout << *dotplot << endl;
  
//  	Alignandum * pro1 = makeProfile( "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA", 1);
//  	Alignandum * pro2 = makeProfile( "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA", 1);
	
//  	pro1->useSegment(10,20);
//  	pro2->useSegment(10,20);
	
//  	// make a Dottor
//  	alignlib::Alignator * dottor = makeAlignatorDummy( dotplot );
	
//  	// make the alignator
//  	alignlib::Alignator * alignator = makeAlignatorDotsSquared( -4, -0.4, dottor );
	
//  	// setup up map_source2dest; m2 will be in col
//  	Alignata *map_source2dest = makeAlignataSet();
	
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
// 	Alignata * ali = makeAlignataSet();
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
//  	Alignata * a = makeAlignataMatrixRow();

//   	Alignandum * seq1 = makeSequence( "XXGXXGXXXXAXXAX" );
//   	Alignandum * seq2 = makeSequence( "XXGXXGXAXXAXXX" );
//  	cout << "seq1: " << *seq1 << endl;
//  	cout << "seq2: " << *seq2 << endl;
//  	Alignator * dottor = makeAlignatorIdentity();
//  	dottor->Align( seq1, seq2, a);
//  	cout << *a << endl;

//  	Alignator * alignator = makeAlignatorDotsDiagonals( 0, -0.1, dottor);
//  	Alignata * ali = makeAlignataSet();
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
