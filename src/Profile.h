/*
  alignlib - a library for aligning protein Profiles

  $Id: ImplProfile.h,v 1.2 2004/01/07 14:35:36 aheger Exp $

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

#ifndef PROFILE_H
#define PROFILE_H 1

#if HAVE_CONFIG_H
#include <config.h>
#endif

#include "alignlib_fwd.h"
#include "Alignandum.h"
#include <iosfwd>

namespace alignlib 
{

    /** A class for profile, that are to be aligned. 
     * 
     * Instances of this class are created by factory functions. This class adds
     * functions that are specific to Profiles to the @ref Alignandum interface.
     * 
     * This class is a protocol class and as such defines only the general interface.
     * 
     * @author Andreas Heger
     * @version $Id$
     * @short A Profile
     */

class Profile : public virtual Alignandum 
{
	
	/* friends ---------------------------------------------------------------------------- */
	friend  std::ostream & operator<<( std::ostream &, const Profile &);

public:
	/** empty constructor */
	Profile();

	/** copy constructor */
	Profile( const Profile &);

	/** destructor */
	virtual ~Profile();

	/** set the @ref Weightor.
	 */
	virtual void setWeightor( const HWeightor & weightor ) = 0;

	/** set the @ref LogOddor.
	 */
	virtual void setLogOddor( const HLogOddor & logoddor ) = 0;

	/** set the @ref Regularizor.
	 */
	virtual void setRegularizor( const HRegularizor & regularizor ) = 0;

	/** get the @ref LogOddor.
	 */
	virtual HWeightor getWeightor() const = 0;
	
	/** get the @ref LogOddor.
	 */
	virtual HLogOddor getLogOddor() const = 0;

	/** get the @ref Regularizor.
	 */
	virtual HRegularizor getRegularizor() const = 0;


	
};

}


#endif /* _Profile_H */

