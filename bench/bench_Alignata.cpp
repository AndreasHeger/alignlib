/*
  alignlib - a library for aligning protein sequences

  $Id: bench_Alignata.cpp,v 1.4 2005/02/24 11:07:25 aheger Exp $

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
#include "Alignata.h"
#include "AlignataIterator.h"
#include "HelpersAlignata.h"

using namespace std;
using namespace alignlib;

#ifdef WITH_DMALLOC
#include <dmalloc.h>
#endif


#define NUM_ITERATIONS 1000
#define NUM_ITERATIONS_LARGE 100000

//-------------------------------------> Initialisation functions <-----------------------------------
void InitCreate( Alignata * a) {
  for (int i = 1; i < 1000; i++) 
    a->addPair(new alignlib::ResiduePAIR( i, i, 1.0));

  a->getLength();
}


//-------------------------------------> Initialisation functions <-----------------------------------
void InitCreateGaps( Alignata * a) {

  for (int i = 1; i < 200; i++) 
    a->addPair(new alignlib::ResiduePAIR( i, i, 1.0));

  for (int i = 300; i < 400; i++) 
    a->addPair(new alignlib::ResiduePAIR( i, i, 1.0));

  for (int i = 500; i < 600; i++) 
    a->addPair(new alignlib::ResiduePAIR( i, i, 1.0));

  for (int i = 800; i < 1000; i++) 
    a->addPair(new alignlib::ResiduePAIR( i, i, 1.0));

  a->getLength();
}

//-------------------------------------> Postprocessing functions <-----------------------------------
void PostClear( Alignata * a) {
  a->clear();
}

//-------------------------------------> Benchmarking functions <-----------------------------------
void BenchmarkAddPair( Alignata * a) {

  for (int i = 1; i < 1000; i++) 
    a->addPair(new  alignlib::ResiduePAIR( i, i, 1.0));
  a->getLength();

}

//-------------------------------------> Benchmarking functions <-----------------------------------
void BenchmarkIterate( Alignata * a) {
  ResiduePAIR p;

  AlignataIterator it(a->begin());
  AlignataIterator it_end(a->end());

  long x;
  for (; it != it_end; it++) {
    p = *it;
    x+=p.mRow;
  }
}

double tval( struct timeval *timval0, struct timeval *timval1) {
  double t;
  
  double t0, t1;
  
  t0 = timval0->tv_sec + timval0->tv_usec*.000001;
  t1 = timval1->tv_sec + timval1->tv_usec*.000001;
  t = t1 - t0;
 
  return t;
}
 
double Benchmark( Alignata * a, 
		  long iterations,
		  void (*ptr_dofunc)(Alignata *), 
		  void (*ptr_prefunc)(Alignata *) = NULL, 
		  void (*ptr_postfunc)(Alignata *) = NULL,
		  void (*ptr_initfunc)(Alignata *) = NULL,
		  void (*ptr_clearfunc)(Alignata *) = NULL) {
  
  struct timeval start_time, finish_time;

  double elapsed = 0;

  if (ptr_initfunc != NULL)
    (*ptr_initfunc)(a);

  for (int x = 0; x < iterations; x++) {
    
    if (ptr_prefunc != NULL)
      (*ptr_prefunc)(a);

    gettimeofday(&start_time,NULL);
  
    (*ptr_dofunc)(a);
  
    gettimeofday(&finish_time,NULL);
    
    elapsed+= tval(&start_time,&finish_time);
    
    if (ptr_postfunc != NULL) 
      (*ptr_postfunc)(a);
  }    

  if (ptr_clearfunc != NULL)
    (*ptr_clearfunc)(a);
 
  return elapsed;
}    

void BenchmarkAll( long iterations,
		   void (*ptr_dofunc)(Alignata *), 
		   void (*ptr_prefunc)(Alignata *) = NULL, 
		   void (*ptr_postfunc)(Alignata *) = NULL,
                   void (*ptr_initfunc)(Alignata *) = NULL,
		   void (*ptr_clearfunc)(Alignata *) = NULL) {

  
  {
    alignlib::Alignata * a1 = alignlib::makeAlignataVector();
    cout << Benchmark( a1, NUM_ITERATIONS, ptr_dofunc, ptr_prefunc, ptr_postfunc, ptr_initfunc, ptr_clearfunc ) << " s\t";
    delete a1;
  }

  {
    alignlib::Alignata * a1 = alignlib::makeAlignataSet();
    cout << Benchmark( a1, NUM_ITERATIONS, ptr_dofunc, ptr_prefunc, ptr_postfunc, ptr_initfunc, ptr_clearfunc ) << " s\t";
    delete a1;
  }

  {
    alignlib::Alignata * a1 = alignlib::makeAlignataHash();
    cout << Benchmark( a1, NUM_ITERATIONS, ptr_dofunc, ptr_prefunc, ptr_postfunc, ptr_initfunc, ptr_clearfunc ) << " s\t"; 
    delete a1;
  }

  {
    alignlib::Alignata * a1 = alignlib::makeAlignataHashDiagonal();
    cout << Benchmark( a1, NUM_ITERATIONS, ptr_dofunc, ptr_prefunc, ptr_postfunc, ptr_initfunc, ptr_clearfunc ) << " s\t";
    delete a1;
  }

  cout << endl;
  
}

int main () {

  cout << "AddPair\t"; BenchmarkAll( NUM_ITERATIONS, &BenchmarkAddPair, NULL, &PostClear, NULL, NULL );
  cout << "Iterate\t"; BenchmarkAll( NUM_ITERATIONS_LARGE, &BenchmarkIterate, NULL, &PostClear, &InitCreate, NULL );
  cout << "IterateGaps\t"; BenchmarkAll( NUM_ITERATIONS_LARGE, &BenchmarkIterate, NULL, &PostClear, &InitCreateGaps, NULL );

//   cout << "---------------------Testing AlignataHash----------------------------------" << endl;
//   a1 = alignlib::makeAlignataHash();
//   TestAlignata( a1 );
//   delete a1;

//   cout << "---------------------Testing AlignataHashDiagonal------------------------------" << endl;
//   a1 = alignlib::makeAlignataHashDiagonal();
//   TestAlignata( a1 );
//   delete a1;

//   cout << "---------------------Testing AlignataMatrixRow-------------------------------" << endl;
//   a1 = alignlib::makeAlignataMatrixRow();
//   TestAlignata( a1 );
//   delete a1;

//   cout << "---------------------Testing AlignataMatrixDiagonal-------------------------------" << endl;
//   a1 = alignlib::makeAlignataMatrixDiagonal();
//   TestAlignata( a1 );
//   delete a1;

//   cout << "---------------------Testing AlignataMatrixUnsorted-------------------------------" << endl;
//   a1 = alignlib::makeAlignataMatrixUnsorted();
//   TestAlignata( a1 );
//   delete a1;

}
