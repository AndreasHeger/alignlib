/*
  alignlib - a library for aligning protein sequences

  $Id: ImplProfile.h,v 1.3 2004/01/07 14:35:35 aheger Exp $

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

#ifndef IMPL_PROFILE_H
#define IMPL_PROFILE_H 1

#include <iosfwd>
#include "alignlib_fwd.h"
#include "alignlib_fwd.h"

#include "ImplAlignandum.h"
#include "Matrix.h"

namespace alignlib 
{

/** 
    Basic implementation for (sequence) profiles, i.e. position specific scoring matrices.
    A profile stores counts, frequencies, and profile-scores for a
    profile and uses different helper objects for the calculation
    of a profile.

    Masking for a profile means:
    1.  profile- value gets set to MASK_VALUE
    2.  counts and frequencies get set to 0

    @author Andreas Heger
    @version $Id: ImplProfile.h,v 1.3 2004/01/07 14:35:35 aheger Exp $
    @short base class for profiles
 */

class ImplProfile : public ImplAlignandum 
{

	friend HAlignandum  extractProfileBinaryCounts( std::istream & input, 
			const Position length,
			const Regularizor *,
			const LogOddor *);

	friend HAlignandum  extractProfileBinaryCountsAsInt( std::istream & input, 
			const Position length,
			int bytes, 
			float scale_factor,
			const Regularizor *,
			const LogOddor *);

	friend HAlignandum addProfile2Profile( HAlignandum dest, 
			const HAlignandum source, 
			const HAlignment map_source2dest );

	friend HAlignandum addSequence2Profile( HAlignandum dest, 
			const HAlignandum source, 
			const HAlignment map_source2dest );

	friend HAlignandum substituteProfileWithProfile( HAlignandum dest, 
			const HAlignandum source, 
			const HAlignment map_source2dest );

	friend HAlignandum rescaleProfileCounts( HAlignandum dest, double scale_factor);

	friend HAlignandum normalizeProfileCounts( HAlignandum dest, Count total_weight );

	friend HAlignandum resetProfile( HAlignandum dest, Position length );

	friend HAlignandum makeProfile( const Position & length,
			const HEncoder &,
			const HRegularizor &,
			const HLogOddor &);

	// friend ProfileFrequencies * exportProfileFrequencies( HAlignandum dest );

public:
	/* constructors and desctructors------------------------------------------------------- */

	/** constructor */
	ImplProfile( 
			const HEncoder & translator,
			const HWeightor & weightor,
			const HRegularizor & regularizor, 
			const HLogOddor & logoddor );	

	ImplProfile( 
			const HMultipleAlignment & src,
			const HEncoder & translator,
			const HWeightor & weightor,
			const HRegularizor & regularizor, 
			const HLogOddor & logoddor );	

	ImplProfile( 
			const Position & length,
			const HEncoder & translator,
			const HWeightor & weightor,
			const HRegularizor & regularizor, 
			const HLogOddor & logoddor );	
	
	/** copy constructor */
	ImplProfile( const ImplProfile &);

	/** destructor */
	virtual ~ImplProfile();

	/** return an identical copy of this object */
	virtual HAlignandum getClone() const;

	/** get internal representation of residue in position pos */
	virtual Residue asResidue( Position pos ) const;

	/* Mutators ------------------------------------------------------------------------------ */

	/** load data into cache, if cacheable type. In this implementation of profiles, this
	will also calculate the profile */
	virtual void prepare() const;						

	/** discard cache, if cacheable type. In this implementation of profiles, this will
     delete the frequencies and profile-types. Only the counts are stored, so that the
     profile can be reconstituted. */
	virtual void release() const;					       

	/** write important member data in a minimally formated way to a stream. Important
	in this sense means data that the user is interested, not internal state variables,
	that would be needed to accurately reconstitute the object.
	Use different "factory" functions to format the output in a way, that you would 
	like to have it (see writeSequenceFasta(...) for an example)
	 */
	virtual void write( std::ostream & output ) const;

	/** swap two positions in the profile */
	virtual void swap( const Position & x, const Position & y );

	/** mask a column */
	virtual void mask( const Position & pos);

	/** save state of object into stream
	 */
	virtual void load( std::istream & input ) ;    

	/** export member data for Scorer and filler objects */
	virtual ScoreMatrix * getScoreMatrix() const ;
	virtual FrequencyMatrix * getFrequencyMatrix() const ;
	virtual CountMatrix * getCountMatrix() const;

protected:

	/** allocate counts for a segment */
	template< class T>
	Matrix<T> * allocateSegment( Matrix<T> * data ) const;

	/** get row with maximum value per column */
	template <class T>
	Residue getMaximumPerColumn( const Matrix<T> * data, const Position & column ) const;

	/** mask a column */
	template<class T>
	void setColumnToValue( Matrix<T> * data, 
			const Position & column,
			const T & value );

	/** write a part of a profile */
	template<class T>
	void writeSegment( std::ostream & output, const Matrix<T> * data ) const;

	/** get residue with most counts in column */
	virtual Residue	getMaximumCount( Position column ) const ;

	/** get residue with highest positive profile score in column */
	virtual Residue	getMaximumScore( Position column ) const ;

	/** get residue with highest frequency in column */
	virtual Residue	getMaximumFrequency( Position column ) const ;

	/** allocate memory for counts 
	 * 
	 * Also sets the width of the profile.
	 * */
	virtual void allocateCounts() const;

	/** allocate memory for frequencies */
	virtual void allocateFrequencies() const;

	/** allocate memory for the profile */
	virtual void allocateScores() const;

	/** fill count matrix */
	virtual void fillCounts( const HMultipleAlignment & src );

	 /** re-set the length of the object
	  * 
	  * This method allocates the memory needed. 
	  */
	 virtual void resize( Position length );
	
	/** save state of object into stream
	 */
	virtual void __save( std::ostream & output, MagicNumberType type = MNNoType ) const;

	/** width of a profile column 
	 * 
	 * This width is set from the alphabet size of the translator and is set if 
	 * allocateCounts() is called.
	 * */
	mutable Residue mProfileWidth;

	/** weightor to use to convert residues to counts 
	 */
	const HWeightor mWeightor;
	
	/** pointer to weighter to use for weighting sequences */
	const HRegularizor mRegularizor;

	/** pointer to objects used for calculating log odds scores */
	const HLogOddor mLogOddor;

	/** pointer to the location of the counts stored in memory */
	mutable CountMatrix *mCountMatrix;			

	/** pointer to the location of the frequencies stored in memory */
	mutable FrequencyMatrix *mFrequencyMatrix;		

	/** pointer to the location of the profile stored in memory */
	mutable ScoreMatrix *mScoreMatrix;		

};

// handle definition for down-casting
typedef boost::shared_ptr<ImplProfile>HImplProfile;



}

#endif /* _PROFILE_H */

