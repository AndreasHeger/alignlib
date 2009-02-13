 /*
  alignlib - a library for aligning protein sequences

  $Id$

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


#include <iostream>
#include "alignlib_fwd.h"
#include "alignlib_interfaces.h"
#include "AlignlibDebug.h"

#include "ImplMultipleAlignatorSimple.h"
#include "HelpersAlignment.h"
#include "HelpersAlignandum.h"
#include "MultAlignmentFormat.h"
using namespace std;

namespace alignlib
{

	HMultipleAlignator makeMultipleAlignatorSimple(
			const HAlignator & alignator )
	{
		debug_func_cerr(5);
		return HMultipleAlignator( new ImplMultipleAlignatorSimple(alignator));
	}


  //----------------------------------------------------------------------------------------
  ImplMultipleAlignatorSimple::ImplMultipleAlignatorSimple(
		  const HAlignator & alignator ) :
		mAlignator( alignator )
  {
	  debug_func_cerr( 5 );
  }

  ImplMultipleAlignatorSimple::~ImplMultipleAlignatorSimple()
  {
      debug_func_cerr(5);
  }

  ImplMultipleAlignatorSimple::ImplMultipleAlignatorSimple( const ImplMultipleAlignatorSimple & src ) :
	  ImplMultipleAlignator(src), mAlignator(src.mAlignator)
  {
  }

  HMultipleAlignator ImplMultipleAlignatorSimple::getClone() const
  {	return HMultipleAlignator( new ImplMultipleAlignatorSimple(*this)); }

  //--------------------------------------------------------------------------
  void ImplMultipleAlignatorSimple::align(
		  HMultAlignment & result,
		  const HAlignandumVector & hsequences ) const
  {
      debug_func_cerr(5);

      AlignandumVector & sequences(*hsequences);

	  result->clear();

	  if (sequences.size() == 0) return;

	  // add the first sequence
	  HAlignment ali(makeAlignmentVector());
	  ali->addDiagonal( 0, sequences[0]->getLength(), 0);
	  result->add( ali );
	  HAlignandumVector aligned(new AlignandumVector());
	  aligned->push_back( sequences[0] );

	  // align the other sequences, expanding the mali
	  // as we go along.
	  for (int x = 1; x < sequences.size(); ++x)
	  {
		  result->expand( aligned );
		  HAlignandum profile( makeProfile( result, aligned ));
		  HAlignment ali(makeAlignmentVector());
		  mAlignator->align( ali, profile, sequences[x]);
		  result->add( ali );
		  aligned->push_back( sequences[x] );
	  }
	  result->expand(aligned);
  }

} // namespace alignlib
