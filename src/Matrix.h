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
#include <cassert>

namespace alignlib 
{

/** A simple fixed-width matrix
 */
template <class T> 
class Matrix 
{

public:
	Matrix (unsigned int r, unsigned int c, T default_value = 0) 
	{
		mRows = r;
		mCols = c;
		mSize = mRows * mCols;
		mMatrix = new T [mSize];
		for (unsigned int i = 0; i < mSize; i++) 
			mMatrix[i] = default_value;
	}                

	T* operator[](unsigned int row) 
	{ 
		assert( row < mRows && row >= 0);
		return &mMatrix[(row * mCols)]; 
	}

	const T* operator[](unsigned int row) const 
	{
		assert( row < mRows );
		return (&mMatrix[(row * mCols)]); 
	}

	T getValue( unsigned int row, unsigned int col) const
	{
		assert( row < mRows);
		assert( col < mCols);
		return mMatrix[(row * mCols) + col ];
	}

	void setValue( unsigned int row, unsigned int col, const T & value)
	{
		assert( row < mRows );
		assert( col < mCols );
		mMatrix[(row * mCols) + col] = value;
	}

	// copy constructor 
	Matrix (const Matrix <T>& src) : 
		mRows(src.mRows), mCols(src.mCols), mSize(src.mSize) 
		{
		mMatrix  = new T [mSize];
		memcpy( mMatrix, src.mMatrix, sizeof(T) * mSize );		
		}

	// assignment operator overloading
	Matrix& operator= (const Matrix <T>&  src) 
	{
		if (src == *this) return (*this) ; // protect against self assignment      
		mRows=src.mRows;
		mCols=src.mCols;
		mSize=mRows*mCols;

		mMatrix  = new T [mSize];
		memcpy( mMatrix, src.mMatrix, sizeof( T) * mSize );

		return (*this);
	}

	// destructor
	~Matrix() { delete [] mMatrix; }

	// equality operator, check element-wise
	bool operator==(const Matrix<T>& other) const 
	{
		if (mRows != other.mRows) return false;
		if (mCols != other.mCols) return false;
		for (unsigned int i = 0; i < mSize; i++)
			if (mMatrix[i] != other.mMatrix[i]) return false;
		return true;
	}

	unsigned int getNumRows() const { return mRows; }
	unsigned int getNumCols() const { return mCols; }

	/** return pointer to location of matrix data 
	 */
	const T * getData()  const
	{
		return mMatrix;
	}

	/** set new data matrix */
	void setData( T * matrix ) 
	{
		mMatrix = matrix;
	}
	
	/** copy data from a data location. This functions copies
	 * as many bytes as it needs. */
	void copyData( T * matrix )
	{
		memcpy( mMatrix, matrix, sizeof(T) * mSize );
	}
	
private:

	/** data location */
	T * mMatrix;
	/** rows in matrix */
	unsigned int mRows;
	/** columns in matrix */
	unsigned int mCols;
	/** size of matrix */
	unsigned int mSize;

};

template<class T> 
std::ostream & operator<< (std::ostream & out, const Matrix<T> & src) 
{

	for (unsigned int i = 0; i < src.getNumRows(); i++) 
	{
		out << i << "=\t";
		for (unsigned int j = 0; j < src.getNumCols(); j++) 
		{
			out << src[i][j] << "\t";
		}
		out << std::endl;
	}
	return out;
}

}
#endif
