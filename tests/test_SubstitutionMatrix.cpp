/*
  alignlib - a library for aligning protein sequences

  $Id: test_HelpersAlignment.cpp,v 1.4 2004/06/02 12:11:38 aheger Exp $

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

using namespace std;

#ifdef WITH_DMALLOC
#include <dmalloc.h>
#endif

using namespace alignlib;

int main ()
{

  Alignment * a = makeAlignmentVector();

  Position from1= 33;
  Alignandum * seq1 = makeSequence("MRDGSLKGLFCTAPKSIRSVRRNSGSFNVTAFTEEQEALVVKSWNAMKKNSGELALKFFLRIFEIAPSAKKLFTFLKDSDIPVEQNPKLKPHATTVFVMTCESAVQLRKAGKVTVRESNLKDLGATHFKYGVADEHFEVTKYALLEPIKEPVPEMWSPELKNAWAEAYDQLAAAIKIEMKPPS");

  Position from2= 13;  
  Alignandum * seq2 = makeSequence("GGTLAIQAQGDLTLAQKKIVRKTWHQLMRNKTSFVTDVFIRIFAYDPSAQNKFPQMAGMSASQLRSSRQMQAHAIRVSSIMSEYVEELDSDILPELLATLARTHDLNKVGADHYNLFAKVLMEALQAELGSDFNEKTRDAWAKAFSVVQAVLLVKHGN");

  std::string ali1 = "+46-1+100";
  std::string ali2 = "+78-4+65";

  fillAlignmentCompressed( a, from1, ali1, from2, ali2);
  rescoreAlignment( a, seq1, seq2 );

  std::cout << *a << std::endl;
  
  a->clear();

  fillAlignmentCompressed( a, from2, ali2, from1, ali1);
  rescoreAlignment( a, seq2, seq1 );

  std::cout << *a << std::endl;  
  
  exit(EXIT_SUCCESS);
  
}
 
