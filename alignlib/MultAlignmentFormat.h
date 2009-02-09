#ifndef Mult_ALIGNMENT_FORMAT_H_
#define Mult_ALIGNMENT_FORMAT_H_

//--------------------------------------------------------------------------------
// Project alignlib
//
// Copyright (C) 2000 Andreas Heger All rights reserved
//
// Author: Andreas Heger <heger@ebi.ac.uk>
//
// $Id$
//--------------------------------------------------------------------------------

#if HAVE_CONFIG_H
#include <config.h>
#endif

#include <iosfwd>
#include <string>

#include "alignlib_fwd.h"

namespace alignlib
{

/**
 *
 * @defgroup MulitpleAlignmentFormats MultAlignment formats.
 * @{
 *
 * MulitpleAlignment formats convert @ref MulitpleAlignment objects into various text based formats and back.
 * Member data is publicly available for easy access in your own formatting functions.
 *
 * Usage example:
 * @code
 *
 * @endcode
 */

/** Base class for Mult alignment formats.
 *
 * 	This is a convenience structure for importing/exporting
 *	Mult alignments.
 *
 *  This class keeps track of alignment coordinates.
 *
 *  @author Andreas Heger
 *  @version $Id$
 *  @short Base class for Mult alignment formats.
 *
 */
struct MultAlignmentFormat
{
	// class member functions
	friend std::ostream & operator<<( std::ostream &, const MultAlignmentFormat &);
	friend std::istream & operator>>( std::istream &, MultAlignmentFormat &);

	// constructors and desctructors
	MultAlignmentFormat ();

	MultAlignmentFormat( const HMultAlignment & src);

	MultAlignmentFormat( std::istream & src);

	MultAlignmentFormat( const std::string & src);

	MultAlignmentFormat( const MultAlignmentFormat &);

	virtual ~MultAlignmentFormat();

	/** fill format from Mult alignment
	 *
	 *	@param src Mult alignment to parse
	 */
	virtual void fill( const HMultAlignment & src);

	/** fill alignment from format
	 *
	 * 	@param dest Alignment
	 */
	virtual void copy( HMultAlignment & dest ) const;

	/** save alignment to stream
	 */
	virtual void save( std::ostream & ) const;

	/** load alignment from stream
	 */
	virtual void load( std::istream &);

	/** string representation of the mali */
	std::string mRepresentation;

};

/**
	Plain Mult alignment format.

	The mali is output in rows.

   	@author Andreas Heger
   	@version $Id$
   	@short Plain Mult alignment format

*/
struct MultAlignmentFormatPlain : public MultAlignmentFormat
{
	// constructors and desctructors
	MultAlignmentFormatPlain ();

	MultAlignmentFormatPlain( const HMultAlignment & src);

	MultAlignmentFormatPlain( std::istream & src);

	MultAlignmentFormatPlain( const std::string & src);

	MultAlignmentFormatPlain (const MultAlignmentFormatPlain &);

	virtual ~MultAlignmentFormatPlain ();

	/** fill blocks from alignment
		@param src Alignment to parse
	 */
	virtual void fill( const HMultAlignment & src);

	/** fill Alignment object with blocks
	 * 	@param dest Alignment
	 */
	virtual void copy( HMultAlignment & dest ) const;

};

/**
	HTML formatted output. Residues are colored according to a palette

   	@author Andreas Heger
   	@version $Id$
   	@short Plain Mult alignment format

*/
struct MultAlignmentFormatHTML : public MultAlignmentFormat
{
	// constructors and desctructors
	MultAlignmentFormatHTML ();

	MultAlignmentFormatHTML( const HMultAlignment & src, const HPalette & palette );

	MultAlignmentFormatHTML( std::istream & src);

	MultAlignmentFormatHTML( const std::string & src);

	MultAlignmentFormatHTML (const MultAlignmentFormatHTML &);

	virtual ~MultAlignmentFormatHTML ();

	/** fill blocks from alignment
		@param src Alignment to parse
	 */
	virtual void fill( const HMultAlignment & src, const HPalette & palette);

	/** fill Alignment object with blocks
	 * 	@param dest Alignment
	 */
	virtual void copy( HMultAlignment & dest ) const;

};


/** @} */

}

#endif
