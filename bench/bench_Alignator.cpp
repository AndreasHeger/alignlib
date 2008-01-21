/*
  alignlib - a library for aligning protein sequences

  $Id: bench_Alignment.cpp,v 1.4 2005/02/24 11:07:25 aheger Exp $

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

/** Test alignata objects
*/

#if HAVE_CONFIG_H
#include <config.h>
#endif

#include <sys/types.h>
#include <sys/socket.h>
#include <sys/time.h>      

#include <iostream>
#include <fstream>

#include <time.h> 

#include "alignlib.h"
#include "Alignment.h"
#include "AlignmentIterator.h"
#include "HelpersAlignment.h"
#include "Alignandum.h"
#include "Alignator.h"
#include "HelpersAlignator.h"
#include "HelpersEncoder.h"
#include "HelpersRegularizor.h"
#include "HelpersLogOddor.h"
#include "HelpersWeightor.h"

using namespace std;
using namespace alignlib;

typedef void(*FunctionType)(HAlignator &, HAlignandum &, HAlignandum &, HAlignment &);
typedef void(*InitFunctionType)(HAlignator &, HAlignandum &, HAlignandum &, HAlignment &);

//-------------------------------------> Initialisation functions <-----------------------------------
void InitSeqSeq( 
			HAlignator & a, 
			HAlignandum  & row, 
			HAlignandum & col,
			HAlignment & result)
{
	result = makeAlignmentVector();
	row = makeSequence( "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA" );
	col = makeSequence( "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA" );
}

void InitSeqProf( 
		 HAlignator & a, 
		 HAlignandum & row, 
		 HAlignandum & col,
		 HAlignment & result)
{
	result = makeAlignmentVector();
	row = makeSequence( "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA" );
	col = makeProfile( "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA", 3,
			getDefaultEncoder(),
			getDefaultWeightor(),
			getDefaultRegularizor(),
			getDefaultLogOddor());
}

void InitProfProf( 
		 HAlignator & a, 
		 HAlignandum & row, 
		 HAlignandum & col,
		 HAlignment & result)
{
	result = makeAlignmentVector();
	row = makeProfile( "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA", 3,
			getDefaultEncoder(),
			getDefaultWeightor(),
			getDefaultRegularizor(),
			getDefaultLogOddor());
			
	col = makeProfile( "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA", 3,
			getDefaultEncoder(),
			getDefaultWeightor(),
			getDefaultRegularizor(),
			getDefaultLogOddor());			
}
						

//-------------------------------------> Initialisation functions <-----------------------------------

//-------------------------------------> Postprocessing functions <-----------------------------------
void ClearAll( 
		HAlignator & a, 
		HAlignandum & row, 
		HAlignandum & col,
		HAlignment & result
)		
{
}


//-------------------------------------> Postprocessing functions <-----------------------------------
void ClearAlignment( 
		 HAlignator & a, 
		 HAlignandum & row, 
		 HAlignandum & col,
		 HAlignment & result
)		
{
  result->clear();
}

//-------------------------------------> Benchmarking functions <-----------------------------------
void BenchmarkAlignment( 
			HAlignator & a, 
			HAlignandum & row, 
			HAlignandum & col,
			HAlignment & result)
{
	a->align( result, row, col );
}

double tval( struct timeval *timval0, struct timeval *timval1) 
{
  double t;
  
  double t0, t1;
  
  t0 = timval0->tv_sec + timval0->tv_usec*.000001;
  t1 = timval1->tv_sec + timval1->tv_usec*.000001;
  t = t1 - t0;
 
  return t;
}
 
double Benchmark( HAlignator & a, 
		  long iterations,
		  FunctionType ptr_dofunc, 
		  FunctionType ptr_prefunc = NULL, 
		  FunctionType ptr_postfunc = NULL,
		  InitFunctionType ptr_initfunc = NULL,
		  FunctionType ptr_clearfunc = NULL) {
  
  struct timeval start_time, finish_time;

  double elapsed = 0;

  HAlignandum row;
  HAlignandum col;
  HAlignment result;
  
  if (ptr_initfunc != NULL)
    (*ptr_initfunc)(a, row, col, result);

  row->prepare();
  col->prepare();
  
  for (int x = 0; x < iterations; x++) 
  {
    
    if (ptr_prefunc != NULL)
      (*ptr_prefunc)(a, row, col, result);

    gettimeofday(&start_time,NULL);
  
    (*ptr_dofunc)(a, row, col, result);
  
    gettimeofday(&finish_time,NULL);
    
    elapsed+= tval(&start_time,&finish_time);
    
    if (ptr_postfunc != NULL) 
      (*ptr_postfunc)(a, row, col, result);
  }    

  if (ptr_clearfunc != NULL)
    (*ptr_clearfunc)(a, row, col, result);
 
  return elapsed;
}    

void BenchmarkAll( 
		long iterations,
		HAlignator & alignator )
{

	// benchmark seq vs seq
    cout << Benchmark( alignator, iterations, &BenchmarkAlignment, NULL, &ClearAlignment, &InitSeqSeq, &ClearAll ) << "\t";
 
    // benchmark seq vs profile
    cout << Benchmark( alignator, iterations, &BenchmarkAlignment, NULL, &ClearAlignment, &InitSeqProf, &ClearAll ) << "\t";
    
    // benchmark prof vs profile
    cout << Benchmark( alignator, iterations, &BenchmarkAlignment, NULL, &ClearAlignment, &InitProfProf, &ClearAll ) << "\t";
  
    cout << endl;
  
}

int main ( int argc, char ** argv) 
{
	
	if (argc != 2)
	{
		std::cerr << "please supply the number of iterations.\n" << std::endl;
		exit(EXIT_FAILURE);
	}
	int num_iterations = atoi( argv[1] );
	
	cout << "alignator\tseqseq\tseqprof\tprofprof" << std::endl;
	{
		HAlignator alignator = makeAlignatorDPFull( ALIGNMENT_LOCAL, -10.0, -2.0 );
		cout << "AlignatorDPFull\t"; BenchmarkAll( num_iterations, alignator );
	}
	exit (EXIT_SUCCESS);
}
