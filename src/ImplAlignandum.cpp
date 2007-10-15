/*
  alignlib - a library for aligning protein sequences

  $Id: ImplAlignandum.cpp,v 1.2 2004/01/07 14:35:33 aheger Exp $

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
#include <iomanip>
#include "alignlib.h"
#include "AlignlibDebug.h"
#include "ImplAlignandum.h"
#include "Translator.h"
#include "HelpersTranslator.h"

#ifdef WITH_DMALLOC
#include <dmalloc.h>
#endif

using namespace std;

namespace alignlib {

//--------------------------------------------------------------------------------------
ImplAlignandum::ImplAlignandum() : 
  mFrom(0), mTo(0), mLength(0), mIsPrepared(false) {
}    

//--------------------------------------------------------------------------------------
ImplAlignandum::~ImplAlignandum () 
{
  debug_func_cerr(5);

}

//--------------------------------------------------------------------------------------
ImplAlignandum::ImplAlignandum(const ImplAlignandum & src) 
{
  debug_func_cerr(5);


    mFrom   = src.mFrom;
    mTo     = src.mTo;
    mLength = src.mLength;
    mIsPrepared = src.mIsPrepared;

}

//--------------------------------------------------------------------------------------
void ImplAlignandum::mask( Position from, Position to) {
  Position j;
  for (j = from; j <= to; j++) 
    mask( j );
}

//--------------------------------------------------------------------------------------
Position ImplAlignandum::getLength() const {
  return (mTo - mFrom + 1);
}

//--------------------------------------------------------------------------------------
void ImplAlignandum::useFullLength() {
  mFrom = 1;
  mTo   = mLength;
}

//--------------------------------------------------------------------------------------
void ImplAlignandum::useSegment(Position from, Position to) {
  mFrom = from;
  mTo   = to;
  if (mTo > mLength)
    mTo = mLength;
}

//--------------------------------------------------------------------------------------
Position ImplAlignandum::getFrom() const {
  return mFrom;
}

//--------------------------------------------------------------------------------------
Position ImplAlignandum::getTo() const {
   return mTo;
}

//--------------------------------------------------------------------------------------
Position ImplAlignandum::getOffset( Position pos) const {
  return (pos + mFrom - 1);
}

//--------------------------------------------------------------------------------------
void ImplAlignandum::setTrueLength( Position length) const {
  mLength = length;
}

//--------------------------------------------------------------------------------------
Position ImplAlignandum::getTrueLength() const {
  return mLength;
}

//--------------------------------------------------------------------------------------
bool ImplAlignandum::isPrepared() const {
    return mIsPrepared;
}

//--------------------------------------------------------------------------------------
void ImplAlignandum::setPrepared( bool flag ) const {
    mIsPrepared = flag;
} 

//--------------------------------------------------------------------------------------
char ImplAlignandum::asChar( Position pos ) const {
  return getDefaultTranslator()->decode( asResidue( pos ));
}


//--------------------------------------------------------------------------------------
// use faster implementations in subclasses, if you prefer
std::string ImplAlignandum::asString() const {
  std::string ret_val("");

  for (Position i = 1; i <= getLength(); i++) 
    ret_val += getDefaultTranslator()->decode( asResidue(i) );
  
  return ret_val;
}
   
  
} // namespace alignlib
