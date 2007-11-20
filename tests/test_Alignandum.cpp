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

#include "HelpersAlignandum.h"
#include "HelpersSequence.h"
#include "HelpersProfile.h"
#include <cassert>

using namespace std;

#ifdef WITH_DMALLOC
#include <dmalloc.h>
#endif

using namespace alignlib;

void testAlignandum( Alignandum * a, const std::string & sample )
{
	assert( a->getFrom() == 0);
	assert( a->getTo() == sample.size() ); 
	assert( sample.size() == a->getLength() );
	for ( int x = 0; x < sample.size(); ++x)
		assert( a->asChar(x) == sample[x]);

	assert( a->asString() == sample );

	// check that useSegment does not interfere with output
	Position from = 1;
	Position to = sample.size() -1;
	a->useSegment( from, to);
	assert( a->getFrom() == from);
	assert( a->getTo() == to );
	for ( int x = 0; x < sample.size(); ++x)
		assert( a->asChar(x) == sample[x]);

	assert( a->asString() == sample );

	// check saving/loading from stream.
	{
		ofstream file("test_Alignandum.tmp", ios::binary);
		a->save( file );
		a->save( file );		
		file.close();
	}

	{
		ifstream file("test_Alignandum.tmp", ios::binary) ;
		
		Alignandum * b = NULL;
		int n = 0;
		while ( b = loadAlignandum( file ) ) 
		{
			assert( a->getFrom() == b->getFrom() );
			assert( a->getLength() == b->getLength() );
			assert( a->getTo() == a->getTo() );
			assert( a->asString() == b->asString() );
			++n; 
		}
		assert( n == 2 );
	}
}


int main () 
{

	{
		std::auto_ptr<Alignandum>a(makeSequence( "ACA") );
		testAlignandum( &*a, "ACA" );
	}

	{
		std::auto_ptr<Alignandum>a(makeSequence( "AAAAACCCCCCCCCCAAAAAAAAAAAAAA") );
		testAlignandum( &*a, "AAAAACCCCCCCCCCAAAAAAAAAAAAAA" );    
	}

	{
		std::auto_ptr<Alignandum>a(makeProfile( "ACDACDACD", 3) );
		testAlignandum( &*a, "ACD" );    
	}

	return EXIT_SUCCESS;

}








