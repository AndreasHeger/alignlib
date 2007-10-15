/*
  alignlib - a library for aligning protein sequences

  $Id: ImplTranslator.cpp,v 1.4 2004/09/16 16:02:38 aheger Exp $

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


#include <string.h>
#include <iostream>
#include "alignlib.h"
#include "AlignlibDebug.h"
#include "ImplTranslator.h"

#ifdef WITH_DMALLOC
#include <dmalloc.h>
#endif

using namespace std;

namespace alignlib {

//---------------------------------------------------------< constructors and destructors >--------------------------------------
ImplTranslator::ImplTranslator () : mEncodingTable(0), mDecodingTable(0) {
}

ImplTranslator::~ImplTranslator () {
}
  
ImplTranslator::ImplTranslator (const ImplTranslator & src ) : 
  mEncodingTable( src.mEncodingTable ), mDecodingTable (src.mDecodingTable ) {
}
  
ImplTranslator::ImplTranslator(  const Residue * encoding_table, const char * decoding_table): 
  mEncodingTable( encoding_table ), mDecodingTable( decoding_table ) {
}

//--------------------------------------------------------------------------------------------------------------------------------
char ImplTranslator::decode( const Residue residue) const { 
  return ( (residue == CODE_GAP) ? CHAR_GAP : mDecodingTable[residue] ); 
}

//--------------------------------------------------------------------------------------------------------------------------------
char * ImplTranslator::decode( const Residue * src, int length ) const 
{
  debug_func_cerr(5);

    char * result = new char[length + 1];
    
    int i;
    for (i = 1; i <= length; i++) 
      result[i-1] = (src[i] == CODE_GAP) ? CODE_GAP : mDecodingTable[src[i]]; 
  
    result[length] = '\0'; 
    return result;
}    

//--------------------------------------------------------------------------------------------------------------------------------
Residue ImplTranslator::encode( const char residue) const 
{ 
  return mEncodingTable[residue]; 
}  

//--------------------------------------------------------------------------------------------------------------------------------
Residue * ImplTranslator::encode( const char* src, int length ) const 
{

    Residue * result = new Residue[length + 1];
    
    int i;
    for (i = 1; i <= length; i++)
      result[i] = mEncodingTable[src[i-1]];
 
    return result;
}

//--------------------------------------------------------------------------------------------------------------------------------
bool ImplTranslator::isValidChar( const char query ) const {
    
  int len = strlen( mDecodingTable );
  
  if ( memchr(mDecodingTable, query, len) || query == CHAR_GAP )
      return true;
  else
      return false;
}    

//--------------------------------------------------------------------------------------------------------------------------------
Residue ImplTranslator::getMaskCode() const {
    return CODE_MASK;
}

//--------------------------------------------------------------------------------------------------------------------------------
Residue ImplTranslator::getGapCode() const {
  return CODE_GAP;
}

//--------------------------------------------------------------------------------------------------------------------------------
char ImplTranslator::getMaskChar() const {
    return CHAR_MASK;
}

//--------------------------------------------------------------------------------------------------------------------------------
char ImplTranslator::getGapChar() const {
  return CHAR_GAP;
}


//--------------------------------------------------------------------------------------------------------------------------------

} // namespace alignlib
