/*
  alignlib - a library for aligning protein sequences

  $Id: test_Alignandum.cpp,v 1.3 2004/06/02 12:11:38 aheger Exp $

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
#include <cassert>
#include <time.h> 

#include "alignlib.h"

#include "Alignatum.h"
#include "HelpersAlignatum.h"
#include "Alignment.h"
#include "HelpersAlignment.h"
#include "Alignandum.h"
#include "HelpersSequence.h"

using namespace std;

#ifdef WITH_DMALLOC
#include <dmalloc.h>
#endif

using namespace alignlib;

void testAlignatum( Alignatum * a, std::string sample )
{
	cout << "testing... with " << *a << "..." ;
	
	Position first_pos = (sample.size() > 0) ? 0 : NO_POS;
	Position last_pos = (sample.size() > 0) ? (int)sample.size() : NO_POS;	
	
	assert( a->getString() == sample );
	assert( a->getFrom() == first_pos);
	assert( a->getTo() == last_pos );
	assert( a->getAlignedLength() == sample.size() );	  
	assert( a->getFullLength() == sample.size() );
	
	// add two terminal gaps on either side
	a->addGaps( 2,2 );
	assert( a->getFrom() == first_pos);
	assert( a->getTo() == last_pos );
	assert( a->getFullLength() == sample.size() );
	assert( a->getAlignedLength() == sample.size() + 4);	
	
	// wrap on alignment - add one gaps into the middle
	std::auto_ptr<Alignment> ali(makeAlignmentVector());
	Position pos = a->getAlignedLength() / 2;
	fillAlignmentIdentity( &*ali, 0, pos);
	fillAlignmentIdentity( &*ali, pos, a->getAlignedLength(), 1);
	
	a->mapOnAlignment( &*ali );
		
	assert( a->getFrom() == first_pos);
	assert( a->getTo() == last_pos );
	assert( a->getFullLength() == sample.size() );
	assert( a->getAlignedLength() == sample.size() + 5);	
	assert( a->getStringReference()[pos] == '-' );
	
	cout << "passed" << endl;	
}

int main () {

  Alignatum * a;
  {	
	  std::cout << "testing empty alignment" << std::endl;
	  a = makeAlignatumFromString( "" );
	  testAlignatum( a, "");
	  delete a;
  }	

  {
	  std::cout << "testing non-empty alignment" << std::endl;
	  a = makeAlignatumFromString("ACDEF");  
	  testAlignatum( a, "ACDEF");	
	  delete a;
  }
  
     {
      std::cout << "testing create from Alignandum" << std::endl;
      std::auto_ptr<Alignandum>s(makeSequence( "ACDEFGHIKL" ));
      a = makeAlignatum( &*s );
      testAlignatum( a, "ACDEFGHIKL");
      delete a;
    }
 
       {
         std::cout << "testing create from Alignandum" << std::endl;
         std::auto_ptr<Alignandum>s(makeSequence( "ACDEFGHIKL" ));
         s->useSegment( 4,6 );
         a = makeAlignatum( &*s );
         testAlignatum( a, "ACDEFGHIKL");
         delete a;
       }

  
  return EXIT_SUCCESS;

}

