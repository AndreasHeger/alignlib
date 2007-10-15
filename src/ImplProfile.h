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
#include "alignlib.h"
#include "ImplAlignandum.h"
#include "HelpersProfile.h"

namespace alignlib {

 
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

class Weightor;
class Regularizor;
class LogOddor;
class MultipleAlignment;
class Alignata;

class ImplProfile : public ImplAlignandum {
    
    friend Alignandum *  extractProfileBinaryCounts( std::istream & input, 
						     const Position length,
						     const Regularizor *,
						     const LogOddor *);
    
    friend Alignandum *  extractProfileBinaryCountsAsInt( std::istream & input, 
							  const Position length,
							  int bytes, 
							  float scale_factor,
							  const Regularizor *,
							  const LogOddor *);

    friend Alignandum * addProfile2Profile( Alignandum * dest, 
					    const Alignandum * source, 
					    const Alignata * map_source2dest );

    friend Alignandum * addSequence2Profile( Alignandum * dest, 
					     const Alignandum * source, 
					     const Alignata * map_source2dest );

    friend Alignandum * substituteProfileWithProfile( Alignandum * dest, 
						      const Alignandum * source, 
						      const Alignata * map_source2dest );

    friend Alignandum * rescaleProfileCounts( Alignandum * dest, double scale_factor);

    friend Alignandum * normalizeProfileCounts( Alignandum * dest, Count total_weight );

    friend Alignandum * fillProfile( Alignandum * dest, 
				     const MultipleAlignment * src, 
				     const Weightor * weightor );

    friend Alignandum * resetProfile( Alignandum * dest, Position length );

    friend Alignandum * makeProfile( Position length,
				     const Regularizor *,
				     const LogOddor *);

    friend ProfileFrequencies * exportProfileFrequencies( Alignandum * dest );
    
 public:
    /* constructors and desctructors------------------------------------------------------- */

    /** empty constructor */
    ImplProfile( const Regularizor * regularizor, 
		 const LogOddor * logoddor );	
    
    /** copy constructor */
    ImplProfile( const ImplProfile &);

    /** destructor */
    virtual ~ImplProfile();

    /** return an identical copy of this object */
    virtual ImplProfile * getClone() const;

    /*------------------------------------------------------------------------------------ */
    /** return a pointer to the member data of this sequence */
    virtual const AlignandumDataProfile & getData() const;

    /** get internal representation of residue in position pos */
    virtual Residue asResidue( Position pos ) const;
    
    /** mask column at position x */
    virtual void mask( Position x);

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

    /** read member data that has been output with the Write subroutine.
     */
    virtual void read( std::istream & input );

    /** shuffle object */
    virtual void shuffle( unsigned int num_iterations = 1,
			  Position window_size = 0);
    
 protected:
    /** get residue with most counts in column */
    virtual Residue	   getMaximumCounts( Position column ) const ;
    
    /** get residue with highest positive profile score in column */
    virtual Residue	   getMaximumProfile ( Position column ) const ;

    /** get residue with highest frequency in column */
    virtual Residue	   getMaximumFrequencies ( Position column ) const ;

    /** allocate memory for counts */
    virtual void allocateCounts();

    /** allocate memory for frequencies */
    virtual void allocateFrequencies() const;
    
    /** allocate memory for the profile */
    virtual void allocateProfile() const;
    
    /*
      virtual void PrintProfile()     const;
      virtual void PrintCounts()      const;
      virtual void PrintFrequencies() const;
    */

    // changed to protected to allow export of data
 protected:
    /** pointer to weighter to use for weighting sequences */
    const Regularizor *  mRegularizor;

    /** pointer to objects used for calculating log odds scores */
    const LogOddor *  mLogOddor;

    /** pointer to the location of the counts stored in memory */
    mutable CountColumn *mCounts;			

    /** pointer to the location of the frequencies stored in memory */
    mutable FrequencyColumn *mFrequencies;		

    /** pointer to the location of the profile stored in memory */
    mutable ProfileColumn *mProfile;		

    mutable AlignandumDataProfile mData;
};
 

}

#endif /* _PROFILE_H */

