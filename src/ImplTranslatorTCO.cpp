/*
  alignlib - a library for aligning protein sequences

  $Id: ImplTranslatorTCO.cpp,v 1.2 2004/01/07 14:35:36 aheger Exp $

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

// Actually this is a bit misleading, there is no new class here, just the data
// and the implementation of a factory function.

#include "ImplTranslator.h"

//--------------------------------------< translation table for blosum matrices >--------------------------------------

namespace alignlib {

    // this string has to be null-terminated, since I use strlen to determine its length
    static char * decoding_table     = ";<=>";

    static Residue encoding_table[131] = { 
            CODE_MASK,													  /* 0 */
	    CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, /* 10 */
	    CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, /* 20 */
	    CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, /* 30 */
	    CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, /* 40 */
	    CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, /* 50 */
	    CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK,         0,         1, /* 60 */
	            2,         3, CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, /* 70 */
	    CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, /* 80 */
	    CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, /* 90 */
	    CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, /* 100 */
	    CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, /* 110 */
	    CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, /* 120 */
	    CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK  /* 130 */
	    };

    // export a default translator object for TCO-sequences
    extern const ImplTranslator DEFAULT_TRANSLATOR_TCO( encoding_table, decoding_table ); 

Translator * makeTranslatorTCO() {
    return new ImplTranslator( encoding_table, decoding_table );
}

} // namespace alignlib

