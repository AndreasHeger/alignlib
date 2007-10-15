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

#include "Alignandum.h"
#include "HelpersSequence.h"
#include "HelpersProfile.h"
#include "Alignator.h"
#include "HelpersAlignator.h"
#include "Alignatum.h"
#include "HelpersAlignatum.h"
#include "Alignata.h"
#include "HelpersAlignata.h"
#include "MultipleAlignment.h"
#include "HelpersMultipleAlignment.h"

#ifdef WITH_DMALLOC
#include <dmalloc.h>
#endif

using namespace std;

const char FILE_SEQ[] = "data/test.seq";
const char FILE_FASTA[] = "data/test.fasta";


int main () {

  {
    // create a multiple alignment 
    alignlib::MultipleAlignment * m1 = alignlib::makeMultipleAlignment();
    m1->Add(alignlib::makeAlignatumFromString("-AAAAA-CCAAA-"));
    m1->Add(alignlib::makeAlignatumFromString("AAAAAA-CCAAAA"));
    m1->Add(alignlib::makeAlignatumFromString("-A-AAA-CCA-A-"));
    m1->Add(alignlib::makeAlignatumFromString("AAAAAAA--AAAA"));
    
    cout << *m1 << endl;
    
    // write segments of the alignment and check numbering
    m1->Write( cout, 0, 13 ); cout << endl;
    m1->Write( cout, 1, 12 ); cout << endl;
    m1->Write( cout, 2, 11 ); cout << endl;
    m1->Write( cout, 3, 10 ); cout << endl;
    
    {
      std::cout << "using sstream" << std::endl;
      std::ostringstream temp;
      std::ostream * p = &temp;
      m1->Write( *p ); *p << ends; std::cout << temp.str() <<endl;
    }
    
    // add two multiple alignments
    alignlib::MultipleAlignment * m2 = alignlib::makeMultipleAlignment();
    
    m2->Add( m1 );
    m2->Add( m1 );
    
    cout << *m2 << endl;
    
    delete m1;
    delete m2;
    
    // add two multiple alignments using an alignment between them:
    
    m1 = alignlib::makeMultipleAlignment();
    m1->Add(alignlib::makeAlignatumFromString("123456789"));
    m1->Add(alignlib::makeAlignatumFromString("123456789"));
    m1->Add(alignlib::makeAlignatumFromString("123456789"));
    
    alignlib::Alignata * ali = alignlib::makeAlignataSet();
    
    ali->addPair( new alignlib::ResiduePAIR( 2,2, 1));
    ali->addPair( new alignlib::ResiduePAIR( 3,3, 1));
    ali->addPair( new alignlib::ResiduePAIR( 5,4, 1));
    ali->addPair( new alignlib::ResiduePAIR( 6,5, 1));
    ali->addPair( new alignlib::ResiduePAIR( 7,7, 1));
    ali->addPair( new alignlib::ResiduePAIR( 8,8, 1));
    
    cout << *ali << endl;
    
    m2 = alignlib::makeMultipleAlignment();
    m2 -> Add( m1 );
    cout << *m2 << endl;
    
    m2->Add( m1, ali, true );
    cout << *m2 << endl;
    
    m2->Add( m1, ali, false);
    cout << *m2 << endl;
    
    delete m1;
    delete m2;
    //-------------------------add new objects using an aligment */
    
    m1 = alignlib::makeMultipleAlignment();
    m1->Add(alignlib::makeAlignatumFromString("123456789"));
    m1->Add(alignlib::makeAlignatumFromString("123456789"));
    m1->Add(alignlib::makeAlignatumFromString("123456789"));
    
    m1->Add(alignlib::makeAlignatumFromString("123456789"), ali, true );
    cout << *m1 << endl;
    
    m1->Add(alignlib::makeAlignatumFromString("123456789"), ali, false );
    cout << *m1 << endl;
    
    delete ali;
    delete m1;
    
  }

  // check mali dots
  std::cout << "## checking MultipleAlignmentDots" << std::endl;
  {

    // create a multiple alignment 
    alignlib::MultipleAlignment * m1 = alignlib::makeMultipleAlignmentDots( false );
    alignlib::MultipleAlignment * m2 = alignlib::makeMultipleAlignmentDots( false );    
    
    alignlib::Alignata * a1 = alignlib::makeAlignataVector();
    alignlib::fillAlignataIdentity( a1, 2, 4, 0 );
    alignlib::fillAlignataIdentity( a1, 7, 8, 0 );

    alignlib::Alignata * a2 = alignlib::makeAlignataVector();
    alignlib::fillAlignataIdentity( a2, 1, 5, +1 );

    alignlib::Alignata * a3 = alignlib::makeAlignataVector();
    alignlib::fillAlignataIdentity( a3, 2, 6, -1 );
    alignlib::fillAlignataIdentity( a3, 7, 7, 0 );
    alignlib::fillAlignataIdentity( a3, 8, 8, +1 );
    
    m1->Add(alignlib::makeAlignatumFromString("ABCDGHJL"), a1);
    m1->Add(alignlib::makeAlignatumFromString(".ABCDEFGH"), a2);
    m1->Add(alignlib::makeAlignatumFromString("BCDEFIJKLM"), a3);

    m2->Add(alignlib::makeAlignatumFromString("ABCDGHJL"), a1);
    m2->Add(alignlib::makeAlignatumFromString(".ABCDEFGH"), a2);
    m2->Add(alignlib::makeAlignatumFromString("BCDEFIJKLM"), a3);
    
    cout << *m1 << endl;
    cout << *m2 << endl;    

    alignlib::Alignata * aa1 = alignlib::makeAlignataVector();
    alignlib::fillAlignataIdentity( aa1, 2, 9, 0 );    

    m1->Add( m2, aa1 );

    cout << *m1 << endl;        
    delete a1;
    delete a2;
    delete a3;
    delete m1;

  }
  
  {

    // create a multiple alignment 
    alignlib::MultipleAlignment * m1 = alignlib::makeMultipleAlignment();
    alignlib::MultipleAlignment * m2 = alignlib::makeMultipleAlignmentDots( false );    
    
    alignlib::Alignata * a1 = alignlib::makeAlignataVector();
    alignlib::fillAlignataIdentity( a1, 1, 8, 0 );

    m1->Add(alignlib::makeAlignatumFromString("ABCDGHIJL"), a1);
    m1->Add(alignlib::makeAlignatumFromString("ABCDGHIJL"), a1);

    m2->Add(alignlib::makeAlignatumFromString("ABCDGHIJL"), a1);
    m2->Add(alignlib::makeAlignatumFromString("ABCDGHIJL"), a1);    


    cout << *m1 << endl;
    cout << *m2 << endl;    

    alignlib::Alignata * aa1 = alignlib::makeAlignataVector();
    alignlib::fillAlignataIdentity( aa1, 2, 5 );
    alignlib::fillAlignataIdentity( aa1, 7, 7, 1 );        

    m1->Add( m2, aa1 );

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
    alignlib::MultipleAlignment * m1 = alignlib::makeMultipleAlignmentDots( true );
    
    alignlib::Alignata * a1 = alignlib::makeAlignataVector();

    alignlib::Alignata * a2 = alignlib::makeAlignataVector();
    alignlib::fillAlignataIdentity( a2, 1, 5, +1 );

    alignlib::Alignata * a3 = alignlib::makeAlignataVector();
    alignlib::fillAlignataIdentity( a3, 2, 6, -1 );
    alignlib::fillAlignataIdentity( a3, 7, 7, 0 );
    alignlib::fillAlignataIdentity( a3, 8, 8, +1 );
    
    m1->Add(alignlib::makeAlignatumFromString("ABCDGHJL") );
    m1->Add(alignlib::makeAlignatumFromString("YABCDEFGH"), a2);
    m1->Add(alignlib::makeAlignatumFromString("BCDEFIJKLM"), a3);

    std::cout << *m1 << endl;
    
    delete m1;
    delete a1;
    delete a2;
    delete a3;

  }
  // check mali dots
  std::cout << "## checking MultipleAlignmentDots options" << std::endl;
  {

    // create a multiple alignment 
    alignlib::MultipleAlignment * m1 = alignlib::makeMultipleAlignmentDots( true );
    
    alignlib::Alignata * a1 = alignlib::makeAlignataVector();
    alignlib::fillAlignataIdentity( a1, 2, 4, 0 );
    alignlib::fillAlignataIdentity( a1, 7, 8, 0 );

    alignlib::Alignata * a2 = alignlib::makeAlignataVector();
    alignlib::fillAlignataIdentity( a2, 1, 5, +1 );

    alignlib::Alignata * a3 = alignlib::makeAlignataVector();
    alignlib::fillAlignataIdentity( a3, 2, 6, -1 );
    alignlib::fillAlignataIdentity( a3, 7, 7, 0 );
    alignlib::fillAlignataIdentity( a3, 8, 8, +1 );
    
    m1->Add(alignlib::makeAlignatumFromString("ABCDGHJL"), a1);
    m1->Add(alignlib::makeAlignatumFromString("YABCDEFGH"), a2);
    m1->Add(alignlib::makeAlignatumFromString("BCDEFIJKLM"), a3);

    std::cout << *m1 << endl;
    
    delete m1;
    delete a1;
    delete a2;
    delete a3;

  }

  // check mali dots
  std::cout << "## checking MultipleAlignmentDots options" << std::endl;
  {

    // create a multiple alignment 
    alignlib::MultipleAlignment * m1 = alignlib::makeMultipleAlignmentDots( true, 0 );
    
    alignlib::Alignata * a1 = alignlib::makeAlignataVector();
    alignlib::fillAlignataIdentity( a1, 2, 4, 0 );
    alignlib::fillAlignataIdentity( a1, 7, 8, 0 );

    alignlib::Alignata * a2 = alignlib::makeAlignataVector();
    alignlib::fillAlignataIdentity( a2, 1, 5, +1 );

    alignlib::Alignata * a3 = alignlib::makeAlignataVector();
    alignlib::fillAlignataIdentity( a3, 2, 6, -1 );
    alignlib::fillAlignataIdentity( a3, 7, 7, 0 );
    alignlib::fillAlignataIdentity( a3, 8, 8, +1 );
    
    m1->Add(alignlib::makeAlignatumFromString("ABCDGHJL"), a1);
    m1->Add(alignlib::makeAlignatumFromString(".ABCDEFGH"), a2);
    m1->Add(alignlib::makeAlignatumFromString("BCDEFIJKLM"), a3);

    cout << *m1 << endl;
    
    delete m1;
    delete a1;
    delete a2;
    delete a3;
  }
  
}








