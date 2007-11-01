# Copyright (C) 2000 Andreas Heger <heger@ebi.ac.uk>
#  
# This file is free software; as a special exception the author gives
# unlimited permission to copy and/or distribute it, with or without 
# modifications, as long as this notice is preserved.
# 
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY, to the extent permitted by law; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

AM_CPPFLAGS = -I$(top_srcdir)/src

check_PROGRAMS = test_Fragmentor test_AlignatorDots test_Alignata  test_Alignator \
		test_MultipleAlignment test_Alignandum test_HelpersAlignata test_Statistics \
		test_SubstitutionMatrix test_Iterator test_IterativeAlignment

TESTS = $(check_PROGRAMS)

LDADD = $(top_srcdir)/src/.libs/libalignlib.a -lm 

test_Fragmentor_SOURCES = test_Fragmentor.cpp 
test_AlignatorDots_SOURCES = test_AlignatorDots.cpp 
test_Alignata_SOURCES = test_Alignata.cpp 
test_Alignator_SOURCES = test_Alignator.cpp 
test_MultipleAlignment_SOURCES = test_MultipleAlignment.cpp 
test_Alignandum_SOURCES = test_Alignandum.cpp 
test_Statistics_SOURCES = test_Statistics.cpp
test_HelpersAlignata_SOURCES = test_HelpersAlignata.cpp
test_SubstitutionMatrix_SOURCES = test_SubstitutionMatrix.cpp
test_Iterator_SOURCES = test_Iterator.cpp
test_IterativeAlignment_SOURCES = test_IterativeAlignment.cpp

