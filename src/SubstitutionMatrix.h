/*
  alignlib - a library for aligning protein sequences

  $Id: SubstitutionMatrix.h,v 1.3 2004/03/19 18:23:41 aheger Exp $

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

#ifndef SUBSTITUTION_MATRIX_H
#define SUBSTITUTION_MATRIX_H 1

#include <iosfwd>

#include "alignlib.h"

namespace alignlib {


/**
   @short base definition of exported substitution matrix data.

   @author Andreas Heger
*/
struct SubstitutionMatrixData {};
 
/**
   @short Interface definition for substitution matrices.
   
   This class is an interface to substitution matrices. It exports a const pointer
   to a memory location to where the matrix stored. The class takes responsibility
   for reading, writing and manipulating matrices. The memory location does belong
   to the matrix it points to!
   
   Matrices are stored row-wise.
   
   A default matrix is supplied as the global libarary symbol DEFAULT_SUBSTMATRIX.

   @author Andreas Heger
   @version $Id: SubstitutionMatrix.h,v 1.3 2004/03/19 18:23:41 aheger Exp $
*/

class SubstitutionMatrix 
{
    friend std::ostream & operator<<( std::ostream &, const SubstitutionMatrix &);
    // class member functions
 public:
    /* constructors and desctructors ----------------------------------------------- */
    /** empty constructor */
    SubstitutionMatrix();	       

    /** copy constructor */
    SubstitutionMatrix  (const SubstitutionMatrix &);
        
    /** desctructor */
    virtual ~SubstitutionMatrix ();

    /* accessors------------------------------------------------------------------------------ */
    /** return pointer to continuos location of matrix in memory */
    virtual const Score * getMatrix() const = 0;

    /** get number of rows of matrix */
    virtual int getNumRows() const = 0;
    
    /** get number of columns of matrix */
    virtual int getNumColums() const= 0;

    /** return the score between characters row and col */
    virtual Score getScore( Residue row, Residue col) const = 0;

    /** returns something, that can be used to access the member data. */
    virtual const SubstitutionMatrixData & getData() const = 0;

    /** write data to stream */
    virtual void write(std::ostream & output) const = 0;
    
    /** clear memory */
    virtual void clear() = 0;

};


}

#endif /* SUBSTITUTION_MATRIX_H */

