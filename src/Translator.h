/*
  alignlib - a library for aligning protein sequences

  $Id: Translator.h,v 1.3 2004/03/19 18:23:41 aheger Exp $

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

#ifndef TRANSLATOR_H
#define TRANSLATOR_H 1

#include "alignlib.h"

namespace alignlib 
{
  
/** @short Interface definition for translators, that translate sequence to internal representation.
    
    Base class for Translators. These objects are responsible for translating a residue or
    a string of residues between real world-representation (e.g. ACVD) and the internal representation
    used in this library.

   @author Andreas Heger
   @version $Id: Translator.h,v 1.3 2004/03/19 18:23:41 aheger Exp $
*/
class Translator 
{
    // class member functions
 public:
    // constructors and destructors

    /** constructor */
    Translator() {};

    /** copy constructor */
    // Translator(const Translator &);

    /** destructor */
    virtual ~Translator () {};

    /** translate a string of residues from internal to real-world representation. A copy
	of the translated (null-terminated) string is returned.
	@param src		pointer to string of residues
	@param length	length of string
    */
    virtual std::string decode( const Residue *src, 
    							int length) const = 0;
    
    /** translate a single residue from internal to real-world representation.
     */
    virtual char decode( const Residue src) const = 0;
 
    /** translate at string of residues from real word presentation to internal representation. A
	copy of the translated string is returned.
	@param src		pointer to string of residues
	@param length	length of string
    */
    virtual Residue * encode( const std::string & src ) const = 0;

    /** translate a single residue from real-world to internal representation.
     */
    virtual Residue encode( const char) const = 0;
 
   /** check, if the supplied character is in the alphabet. */
    virtual bool isValidChar( const char ) const = 0;
    
    /** get code used for a masked character. */
    virtual Residue getMaskCode() const = 0;

    /** get internal code used for a gap. */
    virtual Residue getGapCode()  const = 0;
    
    /** get character used for a masked character. */
    virtual char getMaskChar() const = 0;

    /** get character used for a gap. */
    virtual char getGapChar()  const = 0;

    /** get the size of the alphabet - excluding gap and mask characters */
    virtual int getAlphabetSize() const = 0;
};

}

#endif /* _TRANSLATOR_H */

