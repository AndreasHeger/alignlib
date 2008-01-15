/*
  alignlib - a library for aligning protein sequences

  $Id: ImplAlignatumFasta.h,v 1.2 2004/01/07 14:35:34 aheger Exp $

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

#ifndef IMPL_ALIGNATUM_FASTA_H
#define IMPL_ALIGNATUM_FASTA_H 1

#include <iosfwd>
#include <string>
#include "alignlib_fwd.h"
#include "alignlib_fwd.h"
#include "ImplAlignatum.h"

namespace alignlib 
{

/** Objects of this type represent an aligned object. This class is responsible for
    producing a nice, real-world representation of the aligned object. To this end, 
    this class' provides several methods for modifying an aligned string.

    This alignatum class stores a description line together with the sequence.

    @author Andreas Heger
    @version $Id: ImplAlignatumFasta.h,v 1.2 2004/01/07 14:35:34 aheger Exp $
    @short protocol class for aligned objects
*/

class ImplAlignatumFasta : public ImplAlignatum 
{
 public:
    // constructors and desctructors
    /** when building a new Alignatum sequence by cloning and an alignment, this specifies,
	if the new sequence is in row or column
    */

  /** constructor */
  ImplAlignatumFasta  ();
  
  /** copy constructor */
  ImplAlignatumFasta  (const ImplAlignatumFasta & src);
  
  /** destructor */
  virtual ~ImplAlignatumFasta ();
  
  ImplAlignatumFasta (const std::string & description,
		      const std::string & representation, 
		      Position from = 1, 
		      Position to = 0);

    /** return a copy of the object, mapped if so desired 
	@param ali alignment to be used for the mapping. If empty, return identical copy
	@param skip_gaps true, if gaps in other sequence shall be skipped 
	@param is_in_row true, if alignment is in row
    */
  virtual HAlignatum getClone() const;
    
  /** return an empty copy of this object */
  virtual HAlignatum getNew() const;


  /** write object into stream nicely formatted. Segments are addresed by [from,to)
      @param output stream for result
      @param segment_start beginning of segment 
      @param segment_end end of segment 
  */
  virtual void writeRow( std::ostream & output,
		  const HRenderer & renderer,
		  Position segment_start = 0, 
		  Position segment_end = 0 ) const;

  /** write into stream */
  virtual void  write( std::ostream & output ) const;
  
  /** read from stream */
  virtual void  read( std::istream & input );

 private:
  /** the string representation used for single line output */
  std::string mDescription;
};

}

#endif /* ALIGNATUM_H */

