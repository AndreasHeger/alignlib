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

programs=test_Fragmentor test_AlignatorDots test_Alignment  test_Alignator \
		test_MultipleAlignment test_Alignandum test_HelpersAlignment test_Statistics \
		test_SubstitutionMatrix test_Iterator test_Alignatum \
		test_Tree test_Treetor test_PhyloMatrix test_Distor \
		test_AlignedBlocks test_Weightor \
		test_Translator \
		test_LogOddor

check_PROGRAMS = ${programs}

## test programs get automatically built but not installed
noinst_PROGRAMS = ${programs}

TESTS = $(check_PROGRAMS)

LDADD = $(top_srcdir)/src/.libs/libalignlib.a -lm 

test_Fragmentor_SOURCES = test_Fragmentor.cpp 
test_AlignatorDots_SOURCES = test_AlignatorDots.cpp 
test_Alignment_SOURCES = test_Alignment.cpp 
test_Alignatum_SOURCES = test_Alignatum.cpp
test_Alignator_SOURCES = test_Alignator.cpp 
test_MultipleAlignment_SOURCES = test_MultipleAlignment.cpp 
test_Alignandum_SOURCES = test_Alignandum.cpp 
test_Statistics_SOURCES = test_Statistics.cpp
test_HelpersAlignment_SOURCES = test_HelpersAlignment.cpp
test_AlignedBlocks_SOURCES = test_AlignedBlocks.cpp
test_SubstitutionMatrix_SOURCES = test_SubstitutionMatrix.cpp
test_Iterator_SOURCES = test_Iterator.cpp
test_Treetor_SOURCES = test_Treetor.cpp
test_Tree_SOURCES = test_Tree.cpp
test_PhyloMatrix_SOURCES = test_PhyloMatrix.cpp
test_Distor_SOURCES = test_Distor.cpp
test_Weightor_SOURCES = test_Weightor.cpp
test_Translator_SOURCES = test_Translator.cpp
test_LogOddor_SOURCES = test_LogOddor.cpp
