/*
  alignlib - a library for aligning protein sequences

  $Id: Sequence.cpp,v 1.2 2004/01/07 14:35:37 aheger Exp $

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



#include <string>
#include "Sequence.h" 
#include "AlignException.h"
#include "Translator.h"

#ifdef WITH_DMALLOC
#include <dmalloc.h>
#endif

using namespace std;

namespace alignlib {

//--------------------------------------------------------------------------------------
Sequence::Sequence() : IAlignandum() {
}

//--------------------------------------------------------------------------------------
Sequence::Sequence( const char * src, Translator & translator) : IAlignandum() {
    Position length = strlen( src );

    //!! check for correct translation?
    setTrueLength( length );
    mSequence = translator.encode( src, length );
    setPrepared(true );

}

//--------------------------------------------------------------------------------------
Sequence::Sequence( const Sequence & src ) : IAlignandum( src ) 
{
  debug_func_cerr(5);


  if (mSequence != NULL) delete [] mSequence;
  //!! make exception safe
  mSequence = new Residue[src.getTrueLength()];
  memcpy( mSequence, src.mSequence, src.getTrueLength());
}


//--------------------------------------------------------------------------------------
Sequence::~Sequence() 
{
  debug_func_cerr(5);


  if (mSequence != NULL) 
    delete [] mSequence;
}

//--------------------------------------------------------------------------------------
Residue	Sequence::asResidue(Position n) const { 
  return mSequence[getOffset(n)]; 
}

//--------------------------------------------------------------------------------------
const AlignandumDataSequence * Sequence::getData() const {
    return new AlignandumDataSequence( &mSequence[getOffset(1) - 1]);
}

//--------------------------------------------------------------------------------------
void Sequence::prepare() const {
}

//--------------------------------------------------------------------------------------
void Sequence::release() const {
}
    

} // namespace alignlib
