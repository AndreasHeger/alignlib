# Copyright (C) 2000 Andreas Heger <heger@ebi.ac.uk>
#  
# This file is free software; as a special exception the author gives
# unlimited permission to copy and/or distribute it, with or without 
# modifications, as long as this notice is preserved.
# 
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY, to the extent permitted by law; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.


INCLUDES = -I$(top_srcdir)/src 

AM_CXXFLAGS = $(BOOST_CPPFLAGS)

include_HEADERS = alignlib.h

pkginclude_HEADERS = $(HEADERS_ALIGNLIB) \
			AlignException.h \
			AlignlibDebug.h \
			AlignlibMethods.h \
			Matrix.h \
			$(HEADERS_ALIGNATOR) \
			$(HEADERS_ITERATOR) \
			$(HEADERS_SCORER) \
			$(HEADERS_FRAGMENTOR) \
			$(HEADERS_STATISTICS) \
			$(HEADERS_ALIGNMENT) \
			$(HEADERS_ALIGNANDUM) \
			$(HEADERS_MALI) \
			$(HEADERS_ENCODER) \
			$(HEADERS_SUBSTMATRIX) \
		    $(HEADERS_PHYLOLIB)

HEADERS_ALIGNLIB= alignlib.h alignlib_fwd.h alignlib_default.h \
					alignlib_functions.h alignlib_types.h alignlib_interfaces.h

HEADERS_PHYLOLIB=Tree.h HelpersTree.h \
                        ImplTree.h \
                        DistanceMatrix.h HelpersDistanceMatrix.h \
                        ImplDistanceMatrix.h ImplDistanceMatrixSymmetric.h \
                        Distor.h HelpersDistor.h \
                        ImplDistor.h ImplDistorKimura.h ImplDistorClustal.h ImplDistorDummy.h \
                        Treetor.h HelpersTreetor.h \
                        ImplTreetor.h ImplTreetorDistance.h ImplTreetorDistanceLinkage.h \
                        ImplTreetorDistanceNJ.h


HEADERS_ALIGNANDUM =	Alignandum.h HelpersAlignandum.h \
			Weightor.h HelpersWeightor.h \
			Regularizor.h HelpersRegularizor.h \
			LogOddor.h HelpersLogOddor.h \
			ImplAlignandum.h ImplSequence.h ImplProfile.h \
			ImplWeightor.h ImplWeightorHenikoff.h \
			ImplRegularizor.h ImplRegularizorTatusov.h \
			ImplRegularizorDirichlet.h ImplRegularizorDirichletHash.h \
			ImplRegularizorDirichletInterpolate.h ImplRegularizorDirichletPrecomputed.h \
			ImplLogOddor.h ImplLogOddorUniform.h \
			ImplLogOddorGribskov.h ImplLogOddorBackground.h

HEADERS_ALIGNMENT =	Alignment.h HelpersAlignment.h \
			AlignmentIterator.h AlignmentFormat.h \
			ImplAlignment.h \
			ImplAlignmentVector.h ImplAlignmentSorted.h \
			ImplAlignmentMatrix.h ImplAlignmentMatrixRow.h ImplAlignmentMatrixDiagonal.h \
			ImplAlignmentMatrixUnsorted.h \
			ImplAlignmentBlocks.h
			
HEADERS_ALIGNATOR =	Alignator.h HelpersAlignator.h \
			ImplAlignator.h ImplAlignatorDP.h \
			ImplAlignatorDPFull.h \
			ImplAlignatorIterative.h \
			ImplAlignatorDots.h ImplAlignatorDotsWrap.h \
			ImplAlignatorDotsQuick.h ImplAlignatorDotsDiagonal.h\
			ImplAlignatorIdentity.h ImplAlignatorSimilarity.h \
			ImplAlignatorTuples.h \
			ImplAlignatorPrebuilt.h \
			ImplAlignatorFragments.h \
			ImplAlignatorGroupies.h

HEADERS_ITERATOR=	Iterator2D.h \
			HelpersIterator2D.h \
			ImplIterator2D.h \
			ImplIterator2DFull.h ImplIterator2DBanded.h

HEADERS_SCORER=		Scorer.h HelpersScorer.h \
			ImplScorer.h ImplScorerSequenceSequence.h \
			ImplScorerSequenceProfile.h ImplScorerProfileSequence.h \
			ImplScorerProfileProfile.h

HEADERS_FRAGMENTOR=     Fragmentor.h HelpersFragmentor.h \
			ImplFragmentor.h ImplFragmentorDiagonals.h \
			ImplFragmentorIterative.h ImplFragmentorRepetitive.h

HEADERS_STATISTICS=	Statistics.h

HEADERS_MALI=		MultipleAlignment.h HelpersMultipleAlignment.h \
			Renderer.h HelpersRenderer.h \
			Alignatum.h HelpersAlignatum.h \
			ImplMultipleAlignment.h \
			ImplMultipleAlignmentDots.h \
			ImplRenderer.h ImplRendererPalette.h ImplRendererMView.h ImplRendererColumn.h \
			ImplAlignatum.h  

HEADERS_ENCODER= Encoder.h HelpersEncoder.h \
			ImplEncoder.h

HEADERS_SUBSTMATRIX= HelpersSubstitutionMatrix.h

PARTS = AlignException.cpp \
			AlignlibDebug.cpp \
			$(PARTS_ALIGNATOR) \
			$(PARTS_ITERATOR) \
			$(PARTS_SCORER) \
			$(PARTS_FRAGMENTOR) \
			$(PARTS_STATISTICS) \
			$(PARTS_ALIGNMENT) \
			$(PARTS_ALIGNANDUM) \
			$(PARTS_MALI) \
			$(PARTS_ENCODER) \
			$(PARTS_SUBSTMATRIX) \
			$(PARTS_PHYLOLIB) \
			AlignlibDefault.cpp

PARTS_PHYLOLIB= Tree.cpp HelpersTree.cpp \
                        ImplTree.cpp \
                        DistanceMatrix.cpp HelpersDistanceMatrix.cpp \
                        ImplDistanceMatrix.cpp ImplDistanceMatrixSymmetric.cpp \
                        Distor.cpp HelpersDistor.cpp \
                        ImplDistor.cpp ImplDistorKimura.cpp ImplDistorClustal.cpp ImplDistorDummy.cpp \
                        Treetor.cpp HelpersTreetor.cpp \
                        ImplTreetor.cpp ImplTreetorDistance.cpp ImplTreetorDistanceLinkage.cpp \
                        ImplTreetorDistanceNJ.cpp

PARTS_ALIGNANDUM =   	Alignandum.cpp HelpersAlignandum.cpp HelpersProfile.cpp HelpersSequence.cpp \
			Weightor.cpp HelpersWeightor.cpp \
			Regularizor.cpp HelpersRegularizor.cpp \
			LogOddor.cpp HelpersLogOddor.cpp \
			ImplAlignandum.cpp ImplSequence.cpp ImplProfile.cpp \
			ImplWeightor.cpp ImplWeightorHenikoff.cpp\
			ImplRegularizor.cpp ImplRegularizorTatusov.cpp \
			ImplRegularizorDirichlet.cpp ImplRegularizorDirichletHash.cpp \
			ImplRegularizorDirichletInterpolate.cpp ImplRegularizorDirichletPrecomputed.cpp \
			ImplLogOddor.cpp ImplLogOddorUniform.cpp ImplLogOddorGribskov.cpp \
			ImplLogOddorBackground.cpp ImplLogOddorDirichlet.cpp			

PARTS_ALIGNMENT = Alignment.cpp \
			HelpersAlignment.cpp \
			AlignmentFormat.cpp \
			ImplAlignment.cpp \
			ImplAlignmentVector.cpp \
			ImplAlignmentSorted.cpp \
			ImplAlignmentMatrix.cpp ImplAlignmentMatrixRow.cpp \
			ImplAlignmentMatrixDiagonal.cpp ImplAlignmentMatrixUnsorted.cpp \
			ImplAlignmentBlocks.cpp

PARTS_ALIGNATOR =	Alignator.cpp HelpersAlignator.cpp \
			ImplAlignator.cpp ImplAlignatorDP.cpp \
			ImplAlignatorDPFull.cpp \
			ImplAlignatorIterative.cpp \
			ImplAlignatorDots.cpp ImplAlignatorDotsQuick.cpp ImplAlignatorDotsDiagonal.cpp \
			ImplAlignatorDotsWrap.cpp \
			ImplAlignatorIdentity.cpp ImplAlignatorSimilarity.cpp ImplAlignatorTuples.cpp \
			ImplAlignatorPrebuilt.cpp \
			ImplAlignatorFragments.cpp \
			ImplAlignatorGroupies.cpp

PARTS_ITERATOR=		Iterator2D.cpp HelpersIterator2D.cpp \
			ImplIterator2D.cpp \
			ImplIterator2DFull.cpp ImplIterator2DBanded.cpp

PARTS_SCORER=		Scorer.cpp HelpersScorer.cpp \
			ImplScorer.cpp ImplScorerSequenceSequence.cpp \
			ImplScorerSequenceProfile.cpp ImplScorerProfileSequence.cpp \
			ImplScorerProfileProfile.cpp

PARTS_FRAGMENTOR =	Fragmentor.cpp HelpersFragmentor.cpp \
			ImplFragmentor.cpp ImplFragmentorDiagonals.cpp \
			ImplFragmentorIterative.cpp ImplFragmentorRepetitive.cpp

PARTS_MALI =		MultipleAlignment.cpp HelpersMultipleAlignment.cpp \
			Alignatum.cpp HelpersAlignatum.cpp \
			Renderer.cpp HelpersRenderer.cpp \
			ImplMultipleAlignment.cpp \
			ImplMultipleAlignmentDots.cpp \
			ImplRenderer.cpp ImplRendererPalette.cpp \
			ImplRendererMView.cpp ImplRendererColumn.cpp \
			ImplAlignatum.cpp 

PARTS_STATISTICS =    Statistics.cpp

PARTS_ENCODER=	Encoder.cpp \
			HelpersEncoder.cpp ImplEncoder.cpp 

PARTS_SUBSTMATRIX=HelpersSubstitutionMatrix.cpp

lib_LTLIBRARIES = libalignlib.la  

libalignlib_la_SOURCES = $(PARTS)

LIBOBJS= $(PARTS)


