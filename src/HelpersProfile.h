/*
  alignlib - a library for aligning protein sequences

  $Id: HelpersProfile.h,v 1.4 2004/03/19 18:23:40 aheger Exp $

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

#ifndef HELPERS_PROFILE_H
#define HELPERS_PROFILE_H 1

#include <string>
#include "alignlib.h"
#include "Alignandum.h"

namespace alignlib {

    /** Helper functions for class Alignment:
	
	1. factory functions
	
	2. accessor functions for default objects
	
	3. convenience functions
    */

    
    class MultipleAlignment;
    class Weightor;
    class Regularizor;
    class LogOddor;
    class Alignment;
    class Alignandum; 

    // note: the first column is in 1, not in 0
    /** @short data structure for exporting profile data
     */
    struct AlignandumDataProfile : public AlignandumData {
      CountColumn * mCountsPointer;
      FrequencyColumn * mFrequenciesPointer;
      ProfileColumn * mProfilePointer;
    };

    /* -------------------------------------------------------------------------------------------------------------------- */
    /* 1. factory functions */

    /** create an empty profile. 
     */
    Alignandum * makeProfile( const Regularizor * regularizor = NULL,
			      const LogOddor * logoddor = NULL );

    /** create an empty profile with a given length
     */
    Alignandum * makeProfile( Position length,
			      const Regularizor * regularizor = NULL,
			      const LogOddor * logoddor = NULL );
    
    /** create default profile from a NULL-terminated char-array. 
	@param sequence null terminated string to the concatenated sequences
	@param nsequences number of sequences
    */
    Alignandum * makeProfile( const char * sequences, 
			      int nsequences,
			      const Weightor * weightor = NULL, 
			      const Regularizor * regularizor = NULL,
			      const LogOddor * logoddor = NULL );
    
    /** create default profile from a string. 
	@param sequence null terminated string to the concatenated sequences
	@param nsequences number of sequences
    */
    Alignandum * makeProfile( const std::string & sequences, int nsequences );
    
    /** create default profile from a multiple alignment.
	@param mali multiple alignment
    */

    Alignandum * makeProfile( const MultipleAlignment * mali, 
			      const Weightor * weightor = NULL, 
			      const Regularizor * regularizor = NULL,
			      const LogOddor * logoddor = NULL );

    /** read counts in binary format from stream and return a profile */
    Alignandum * extractProfileBinaryCounts( std::istream & input, 
					     const Position max_length = MAX_SEQLEN,
					     const Regularizor * regularizor = NULL,
					     const LogOddor * logoddor = NULL );

    Alignandum * extractProfileBinaryCountsAsInt( std::istream & input, 
						  const Position max_length = MAX_SEQLEN,
						  int bytes = 2, 
						  float scale_factor = 1,
						  const Regularizor * regularizor = NULL,
						  const LogOddor * logoddor = NULL );

    /* -------------------------------------------------------------------------------------------- */
    /* 2. accessor functions for default objects */
    

    /* -------------------------------------------------------------------------------------------- */
    /* 3. convenience functions */

    /** write Profile counts to stream in compressed format. Counts are stored columnwise. */
    void writeProfileBinaryCounts( std::ostream & output, const Alignandum * src);

    /** write Profile counts to stream in compressed format. Counts are stored columnwise as integers. */
    void writeProfileBinaryCountsAsInt( std::ostream & output, 
					const Alignandum * src, 
					int bytes = 2, 
					float scale_factor = 1);

    /** add counts of profile source to profile dest, using the mapping provided, where dest is in col and
	source is in row 
    */
    Alignandum * addProfile2Profile( Alignandum * dest, 
				     const Alignandum * source, 
				     const Alignment * map_source2dest );

    /** add counts of profile source to profile dest, using the mapping provided, where dest is in col and
	source is in row 
    */
    Alignandum * addSequence2Profile( Alignandum * dest, 
				      const Alignandum * source, 
				      const Alignment * map_source2dest );
    
    /** substitutes columns in profile dest by columns in profile row using the mapping provided, 
	where dest is in col and source is in row
     */
    Alignandum * substituteProfileWithProfile( Alignandum * dest, 
					       const Alignandum * source, 
					       const Alignment * map_source2dest );
    
    /** rescale counts from a profile by multiplying each entry by the scale_factor */
    Alignandum * rescaleProfileCounts( Alignandum * dest, double scale_factor );

    /** normalize counts from a profile so that all sum to total_weight per column*/
    Alignandum * normalizeProfileCounts( Alignandum * dest, Count total_weight );

    /** fill a profile with a multiple alignment using a specified weightor */
    Alignandum * fillProfile( Alignandum * dest, 
			      const MultipleAlignment * src, 
			      const Weightor * weightor = NULL );

    /** reset a profile to a new length. Clear old values.
     */
    Alignandum * resetProfile( Alignandum * dest, Position new_length );

    /** export frequencies of a profile. Note, that you receive
	a copy and thus need to free it.
    */
    ProfileFrequencies * exportProfileFrequencies( Alignandum * dest );

}

#endif /* PROFILE_H */

