/*
  alignlib - a library for aligning protein sequences

  $Id: ImplAlignatorFragmentsSquared.cpp,v 1.2 2004/01/07 14:35:34 aheger Exp $

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


#include <iostream>
#include <iomanip>
#include <stdio.h>

#include <map>
#include <vector>
#include <algorithm>

#include "alignlib.h"
#include "alignlib_fwd.h"
#include "AlignlibDebug.h"
#include "AlignException.h"
#include "ImplAlignatorFragmentsSquared.h"
#include "Alignandum.h"

#include "Alignment.h"
#include "HelpersAlignment.h"

#include "HelpersSubstitutionMatrix.h"

using namespace std;

namespace alignlib 
{

#define NOFRAGMENT -1

struct BorderPosition 
{
	BorderPosition() : mRow(0), mTop(false), mFragment(NOFRAGMENT) {}
	BorderPosition( Position row, bool top, Fragment fragment) : 
		mRow(row), mTop(top), mFragment(fragment) {}
	Position mRow;
	bool mTop;
	Fragment mFragment;
};

class BorderComparator 
{
public:
	int operator() (const BorderPosition & x, const BorderPosition & y) const 
	{
		return x.mRow < y.mRow;
	}
};


/*---------------------factory functions ---------------------------------- */

/** make an alignator object, which does a dot-alignment. The default version can be given an AlignmentMatrix-
	object */
HAlignator makeAlignatorFragmentsSquared(
		Score gop, 
		Score gep, 
		HFragmentor & fragmentor) 
{
	return HAlignator( new ImplAlignatorFragmentsSquared( fragmentor, gop, gep, gop, gep));
}

//----------------------------------------------------------------------------------------------------------------------------------------
/** constructors and destructors */
ImplAlignatorFragmentsSquared::ImplAlignatorFragmentsSquared(
		HFragmentor & fragmentor,
		Score row_gop, Score row_gep, 
		Score col_gop, Score col_gep ):
			ImplAlignatorFragments( fragmentor, row_gop, row_gep, col_gop, col_gep ) 
			{
			}

		//------------------------------------------------------------------------------------------
		ImplAlignatorFragmentsSquared::ImplAlignatorFragmentsSquared( const ImplAlignatorFragmentsSquared & src ) : 
			ImplAlignatorFragments( src ) 
			{
			debug_func_cerr(5);

			}

		//------------------------------------------------------------------------------------------
		ImplAlignatorFragmentsSquared::~ImplAlignatorFragmentsSquared() 
		{
			debug_func_cerr(5);

		}

		//----------------------------------------------------------------------------------------------------------
		HAlignator ImplAlignatorFragmentsSquared::getClone() const 
		{
			return HAlignator( new ImplAlignatorFragmentsSquared( *this ) );
		}

		//------------------------------------------------------------------------------------------
		//------------------------------------------------------------------------------------------
		Score ImplAlignatorFragmentsSquared::getGapCost( Dot x1, Dot x2 ) const 
		{

			Position r1 = (*mFragments)[x1]->getRowTo();
			Position r2 = (*mFragments)[x2]->getRowFrom();  
			Position c1 = (*mFragments)[x1]->getColTo();
			Position c2 = (*mFragments)[x2]->getColFrom();  

			Score gap_cost = 0;
			Position d;

			if ((d = (r2 - r1)) > 1)
				gap_cost += mRowGop + d * mRowGep;

			if ((d = (c2 - c1)) > 1)
				gap_cost += mColGop + d * mColGep;

			return gap_cost;
		}

		//-------------------------------------------< Alignment subroutine >----------------------------------------------
		void ImplAlignatorFragmentsSquared::performAlignment( 
				HAlignment & ali,
				const HAlignandum & prow, 
				const HAlignandum & pcol ) 
		{

			/**
     Overview over the algorithm:

     - proceed row-wise.
     - update score of path up to current fragment,
     if left border is touched.
     - add fragment to search-region, if bottom is
     higher than the current row.

			 */
			// #define DEBUG
			// #define DEBUG2

			Fragment global_best_fragment = NOFRAGMENT;
			Score global_best_score = 0;

			// fill borders with fragments
			vector<BorderPosition> borders;
			typedef vector<BorderPosition>::const_iterator BorderIterator;
			for (Fragment i = 0; i < mNFragments; i++) {
				borders.push_back(BorderPosition((*mFragments)[i]->getRowFrom(), true, i));
				borders.push_back(BorderPosition((*mFragments)[i]->getRowTo(), false, i));    
			}

			// sort borders according to row
			sort( borders.begin(), borders.end(), BorderComparator());

			// create data structure for search region
			// sort fragments in search region by increasing column
			typedef multimap <Position, Fragment> MyFragmentSet;
			MyFragmentSet search_region;

			// array with scores of fragments
			vector<Score> scores(mNFragments,0);

			// array with fragments in current row
			vector<Fragment> fragment_stack;

			Position last_row = 0;

			//----------------------------------> main alignment loop <----------------------------------------------------
			BorderIterator it(borders.begin()), it_end(borders.end());
			for ( ; it!=it_end; ++it ) {	   

				Fragment current_fragment = it->mFragment;
				Position current_row = it->mRow;                         
				bool is_top      = it->mTop;

#ifdef DEBUG
				cout << "current_fragment=" << current_fragment
				<< " current_row=" << current_row
				<< " is_top=" << is_top 
				<< endl;
#endif

				/* if a new row is entered, enter fragments from stack to search-area */
				if (current_row != last_row) {

					while (!fragment_stack.empty()) {
						Fragment fragment = fragment_stack.back();
						fragment_stack.pop_back();
						search_region.insert(pair<Position, Fragment>((*mFragments)[fragment]->getColTo(), fragment));
					}
					last_row = current_row;
				}

				if (is_top) {
					/* search search-area: always lookup starting at col 1 until current_col - 1. Try to find
	 a positive trace leading to current fragment. If it were negative, it would not be part of
	 the optimum alignment up to current_fragment. */
					Fragment search_best_fragment   = NOFRAGMENT;
					Position current_col = (*mFragments)[current_fragment]->getColFrom();
					Score search_best_score = 0;

#ifdef DEBUG2
					cout << "SEARCH_AREA" << endl;
					for (MyFragmentSet::iterator search_it = search_region.begin(); search_it != search_region.end(); ++search_it)
						cout << "  [" << (*search_it).first << ", " << (*search_it).second << "]" 
						<< "=" << scores[(*search_it).second] 
						                 << " " << mTrace[(*search_it).second] 
						                                  << endl;
#endif

					MyFragmentSet::const_iterator search_it(search_region.begin()), search_it_end(search_region.end());
					while (search_it != search_it_end && ((*search_it).first < current_col)) 
					{

						Fragment const search_fragment = (*search_it).second;
						Score search_score = scores[search_fragment];

						if (search_score > 0) 
						{
							search_score += getGapCost( search_fragment, current_fragment);
							if (search_score >= search_best_score) 
							{
								search_best_score = search_score;
								search_best_fragment   = search_fragment;
							}
						}

						++search_it;
					}

					/* no positive trace found, new trace starts at current fragment */
					if (search_best_fragment == NOFRAGMENT)
						search_best_score = (*mFragments)[current_fragment]->getScore();
					else
						search_best_score += (*mFragments)[current_fragment]->getScore();

#ifdef DEBUG
					cout << "current_fragment=" << current_fragment 
					<< " current_row=" << current_row 
					<< " current_col=" << current_col 
					<< endl;
					cout << "search_best_fragment=" << search_best_fragment 
					<< " search_best_score=" << search_best_score 
					<< endl;
#endif

					/* do local alignment, traces with score <= 0 are skipped */
					if (search_best_score < 0)
						continue;

					scores[current_fragment] = search_best_score;
					mTrace[current_fragment] = search_best_fragment;

					/* remember end point of best trace */
					if (search_best_score > global_best_score) {
						global_best_score = search_best_score;
						global_best_fragment   = current_fragment;
					}

				} else {
					/* add fragment to list of fragments, that are to be added to the
	 search area, once the current row has finished */
					fragment_stack.push_back(current_fragment);    
				}
			} 
			/* ---->end of alignment loop<----- */

			mLastFragment= global_best_fragment;
			mScore  = global_best_score;

#ifdef DEBUG
			cout << "global_best_fragment=" << global_best_fragment << " global_best_score=" << global_best_score << endl;
#endif

			//--------------> cleaning up <---------------------------------------------------------------
			//--------------------------------------------------------------------------------------------------

		}

} // namespace alignlib






