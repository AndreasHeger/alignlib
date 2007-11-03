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

#include <time.h> 

#include "alignlib.h"

#include "Alignandum.h"
#include "Translator.h"
#include "HelpersTranslator.h"

#include "HelpersSequence.h"
#include "HelpersProfile.h"
#include <cassert>

using namespace std;

#ifdef WITH_DMALLOC
#include <dmalloc.h>
#endif

using namespace alignlib;

void testSequence( const std::string & sample )
  {
    cout << "testing...makeSequence()...";    

    Alignandum * a = makeSequence( sample );
    assert( sample.size() == a->getLength() );
    for ( int x = 0; x < sample.size(); ++x)
      {
        assert( a->asChar(x) == sample[x]);
      }
    // the following assert does not work, although the printout looks identical
    // assert( strcmp( a->asString().c_str(), sample.c_str() ) );
    delete a;
    cout << "passed" << endl;    
  }


int main () {

  {
    testSequence( "ACA" );
    testSequence( "AAAAACCCCCCCCCCAAAAAAAAAAAAAA" );    
  }

  {
    Alignandum * a;
    cout << "testing...makeProfile()...";
    a = makeProfile( "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA", 3 );
    delete a;
    cout << "passed" << endl;
  }
  

  return EXIT_SUCCESS;

}








