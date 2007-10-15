# Copyright (C) 2000 Andreas Heger <heger@ebi.ac.uk>
#  
# This file is free software; as a special exception the author gives
# unlimited permission to copy and/or distribute it, with or without 
# modifications, as long as this notice is preserved.
# 
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY, to the extent permitted by law; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

# check_PROGRAMS = test_objects test_profile test_ma test_alignatum test_tree test_matrix test_alignment
noinst_PROGRAMS = bench_Alignata #   test_Alignator  # test_AlignatorDots # test_Alignata test_MultipleAlignment

LDADD = $(top_srcdir)/src/.libs/libalignlib.a -lm # -ldmallocxx

# test_AlignatorDots_SOURCES = test_AlignatorDots.cpp # dmallocc.cc
# test_Alignator_SOURCES = test_Alignator.cpp # dmallocc.cc
# test_MultipleAlignment_SOURCES = test_MultipleAlignment.cpp # dmallocc.cc
# test_Profile_SOURCES = test_Profile.cpp # dmallocc.cc
bench_Alignata_SOURCES = bench_Alignata.cpp # dmallocc.cc

INCLUDES = -I$(top_srcdir)/src 

AUTOMAKE_OPTIONS = 
