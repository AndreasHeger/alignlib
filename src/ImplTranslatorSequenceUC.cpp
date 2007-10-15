/*
  alignlib - a library for aligning protein sequences

  $Id: ImplTranslatorSequenceUC.cpp,v 1.2 2004/01/07 14:35:36 aheger Exp $

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

/* 
   Translate only upper case characters. Do not use lower case characters. They are treated
   as masked residues.

*/

#include "HelpersTranslator.h"
#include "ImplTranslator.h"

//--------------------------------------< translation table for blosum matrices >--------------------------------------

namespace alignlib {

    // this string has to be null-terminated, since I use strlen to determine its length
    static char * decoding_table = "ACDEFGHIKLMNPQRSTVWYX";	

    // translate upper and lower-case characters
    static Residue encoding_table[131] = { 
            CODE_MASK,													/* 0 */
	    CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, /* 10 */
	    CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, /* 20 */
            CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, /* 30 */
	    CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, /* 40 */
	    CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, /* 50 */
	    CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, /* 60 */
	    CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK,         0, CODE_MASK,         1,         2,         3,         4, /* 70 */
	            5,         6,         7, CODE_MASK,         8,         9,        10,        11, CODE_MASK,        12, /* 80 */
	           13,        14,        15,        16, CODE_MASK,        17,        18, CODE_MASK,        19, CODE_MASK, /* 90 */
	    CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, /* 100 */
	    CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, /* 110 */
	    CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, /* 120 */
	    CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, CODE_MASK, /* 130 */
		};

    
    Translator * makeTranslatorSequenceUC() { return new ImplTranslator( encoding_table, decoding_table ); }
    
} // namespace alignlib

