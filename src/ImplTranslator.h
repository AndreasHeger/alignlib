/*
  alignlib - a library for aligning protein sequences

  $Id: ImplTranslator.h,v 1.3 2004/03/19 18:23:41 aheger Exp $

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

#if HAVE_CONFIG_H
#include <config.h>
#endif

#ifndef IMPL_TRANSLATOR_H
#define IMPL_TRANSLATOR_H 1

#include <string>
#include "Translator.h"

namespace alignlib 
{
    
#define CODE_MASK       20		/* character corresponding to a masked residues, do not change, this position is fixed in some translation tables */
#define CODE_GAP       100		/* code corresponding to gap */

#define UNKNOWN_CODE   'X'		/* character corresponding to a residue with unknown code */
#define CHAR_GAP       '-'		/* character used for a gap */
#define CHAR_MASK      'X'		/* character used for a gap */

/** @short Basic and complete implementation of translators.

    This implementation uses two translating tables, that have to be
    supplied to the constructor.
    
    @author Andreas Heger
    @version $Id: ImplTranslator.h,v 1.3 2004/03/19 18:23:41 aheger Exp $
*/
class ImplTranslator : public Translator 
{
    // class member functions
 public:
    /** constructor */
    ImplTranslator();

    /** create a translator with two pointers to the translation tables 
	@param encoding_table	pointer to array of Residue used for encoding
	@param decoding_table	pointer to array of chars used for decoding
    */
    ImplTranslator( const std::string & alphabet );

    /** copy constructor */
    ImplTranslator(const ImplTranslator &);

    /** destructor */
    virtual ~ImplTranslator ();

    /** translate a string of residues from internal to real world representation
	@param src		pointer to string of residues
	@param length	length of string
    */
    virtual std::string decode( const Residue *src, int length) const;
    
    /** translate a single residue from internal to real world representation.
     */
    virtual char decode( const Residue src) const;
 
    /** translate at string of residues from real word presentation to internal representation.
	@param src		pointer to string of residues
	@param length	length of string
    */
    virtual Residue * encode( const std::string & src) const;

    /** translate a single residue from real world to internal representation.
     */
    virtual Residue encode( const char) const;
 
   /** check, if the supplied character is in the alphabet. */
    virtual bool isValidChar( const char ) const;

    /** get code used for a masked character. */
    virtual Residue getMaskCode() const; 

    /** get internal code used for a gap. */
    virtual Residue getGapCode()  const; 

    /** get character used for a masked character. */
    virtual char getMaskChar() const;

    /** get character used for a gap. */
    virtual char getGapChar()  const;
    
    /** get the size of the alphabet - excluding gap and mask characters */
    virtual int getAlphabetSize() const;
        
 private:
	
	/** the alphabet */
	std::string mAlphabet;
	
    /** pointer to encoding table */
    Residue * mEncodingTable;

};




}

#endif /* IMPL_TRANSLATOR_H */

