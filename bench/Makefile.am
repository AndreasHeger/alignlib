# Copyright (C) 2000 Andreas Heger <heger@ebi.ac.uk>
#  
# This file is free software; as a special exception the author gives
# unlimited permission to copy and/or distribute it, with or without 
# modifications, as long as this notice is preserved.
# 
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY, to the extent permitted by law; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

programs = bench_Alignment bench_Alignator
noinst_PROGRAMS =$(programs) 

LDADD = $(top_srcdir)/alignlib/.libs/libalignlib.a -lm 

bench_Alignment_SOURCES = bench_Alignment.cpp 
bench_Alignator_SOURCES = bench_Alignator.cpp

INCLUDES = -I$(top_srcdir)/alignlib 


