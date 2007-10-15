/*
  alignlib - a library for aligning protein sequences

  $Id: Matrix.h,v 1.2 2004/01/07 14:35:37 aheger Exp $

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

#ifndef MATRIX_H
#define MATRIX_H 1

#include <iostream>

namespace alignlib {

template <class T> class Matrix {

 public:
  T * mMatrix ;
    
 private:
  unsigned int mRows;
  unsigned int mCols;
  unsigned int mSize;
	
 public:
       Matrix (unsigned int r, unsigned int c, T default_value = 0) {
	    mRows = r;
	    mCols =c;
	    mSize = mRows * mCols;
	    mMatrix = new T [mSize];
	    for (unsigned int i = 0; i < mSize; i++) 
		mMatrix[i] = default_value;
	}                

	T* operator[](int x) { return (&mMatrix[(x * mCols)]); }
  
	const T* operator[](int x) const { return (&mMatrix[(x * mCols)]); }
  
	// copy constructor 
	Matrix (const Matrix <T>& from) : mRows(from.mRows), mCols(from.mCols), mSize(from.mSize) {
	    mMatrix  = new T [mSize];
	    for (unsigned int i = 0; i < mSize; i++) 
		mMatrix [i] = from.mMatrix[i];
	}

	// assignment operator overloading
	Matrix& operator= (const Matrix <T>&  from) {
	    if (from == *this) return (*this) ; // protect against self assignment      
	    mRows=from.mRows;
	    mCols=from.mCols;
	    mSize=mRows*mCols;

	    mMatrix  = new T [mSize];
	    for (unsigned int i = 0; i < mSize; i++)
		mMatrix [i] = from.mMatrix[i];

	    return (*this);
	}
  
	// destructor
	~Matrix() { delete [] mMatrix; }

	// equality operator, check element-wise
	bool operator==(const Matrix<T>& other) const {
	  if (mRows != other.mRows) return false;
	  if (mCols != other.mCols) return false;
	  for (unsigned int i = 0; i < mSize; i++)
	    if (mMatrix[i] != other.mMatrix[i]) return false;
	  return true;
	}
	
	unsigned int getNumRows() const { return mRows; }
	unsigned int getNumCols() const { return mCols; }

};

template<class T> 
std::ostream & operator<< (std::ostream & out, const Matrix<T> & src) {
  
 for (unsigned int i = 0; i < src.getNumRows(); i++) {
     out << i << "=\t";
    for (unsigned int j = 0; j < src.getNumCols(); j++) {
	out << src[i][j] << "\t";
    }
    out << std::endl;
  }
 return out;
}

}
#endif
