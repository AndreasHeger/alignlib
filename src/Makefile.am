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

include_HEADERS = alignlib.h

pkginclude_HEADERS = 	alignlib.h \
			AlignException.h \
			AlignlibDebug.h \
			AlignlibMethods.h \
			Matrix.h \
			$(HEADERS_ALIGNATOR) \
			$(HEADERS_ITERATOR) \
			$(HEADERS_SCORER) \
			$(HEADERS_FRAGMENTOR) \
			$(HEADERS_STATISTICS) \
			$(HEADERS_ALIGNATA) \
			$(HEADERS_ALIGNANDUM) \
			$(HEADERS_MALI) \
			$(HEADERS_TRANSLATOR) \
			$(HEADERS_SUBSTMATRIX) \
		    $(HEADERS_PHYLOLIB)

HEADERS_PHYLOLIB=Tree.h HelpersTree.h \
                        ImplTree.h \
                        PhyloMatrix.h HelpersPhyloMatrix.h \
                        ImplPhyloMatrix.h ImplPhyloMatrixSymmetric.h \
                        Distor.h HelpersDistor.h \
                        ImplDistor.h ImplDistorKimura.h ImplDistorClustal.h ImplDistorDummy.h \
                        Treetor.h HelpersTreetor.h \
                        ImplTreetor.h ImplTreetorDistance.h ImplTreetorDistanceLinkage.h \
                        ImplTreetorDistanceNJ.h
                        


HEADERS_ALIGNANDUM =	Alignandum.h HelpersAlignandum.h HelpersSequence.h HelpersProfile.h \
			Weightor.h HelpersWeightor.h \
			Regularizor.h HelpersRegularizor.h \
			LogOddor.h HelpersLogOddor.h \
			ImplAlignandum.h ImplSequence.h ImplProfile.h \
			ImplNoWeightor.h ImplWeightor.h ImplWeightorHenikoff.h ImplWeightorHenikoffKimmen.h \
			ImplRegularizor.h ImplNoRegularizor.h \
			ImplRegularizorDirichlet.h ImplRegularizorDirichletHash.h \
			ImplRegularizorDirichletInterpolate.h ImplRegularizorDirichletPrecomputed.h \
	                ImplLogOddor.h 

HEADERS_ALIGNATA =	Alignata.h HelpersAlignata.h \
			AlignataIterator.h \
			ImplAlignata.h \
			ImplAlignataVector.h \
			ImplAlignataSet.h ImplAlignataSetCol.h \
			ImplAlignataHash.h ImplAlignataHashDiagonal.h \
			ImplAlignataMatrix.h ImplAlignataMatrixRow.h ImplAlignataMatrixDiagonal.h \
			ImplAlignataMatrixUnsorted.h 

HEADERS_ALIGNATOR =	Alignator.h HelpersAlignator.h \
			ImplAlignator.h ImplAlignatorDP.h \
			ImplAlignatorDPFull.h \
			ImplAlignatorDots.h ImplAlignatorDotsWrap.h \
			ImplAlignatorDotsSquared.h ImplAlignatorDotsSquaredDiagonal.h\
			ImplAlignatorIdentity.h ImplAlignatorSimilarity.h \
			ImplAlignatorTuples.h ImplAlignatorDummy.h \
			ImplAlignatorFragments.h ImplAlignatorFragmentsSquared.h \
			ImplAlignatorPublishAlignata.h 

HEADERS_ITERATOR=	Iterator2D.h Iterator2DFull.h \
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
			ImplRenderer.h ImplRendererMView.h ImplRendererColumn.h \
			ImplAlignatum.h ImplAlignatumFasta.h 

HEADERS_TRANSLATOR= 	Translator.h HelpersTranslator.h \
			ImplTranslator.h

HEADERS_SUBSTMATRIX=    SubstitutionMatrix.h HelpersSubstitutionMatrix.h \
			ImplSubstitutionMatrix.h ImplSubstitutionMatrixAA.h

PARTS =               AlignException.cpp \
			AlignlibDebug.cpp \
			$(PARTS_ALIGNATOR) \
			$(PARTS_ITERATOR) \
			$(PARTS_SCORER) \
			$(PARTS_FRAGMENTOR) \
			$(PARTS_STATISTICS) \
			$(PARTS_ALIGNATA) \
			$(PARTS_ALIGNANDUM) \
			$(PARTS_MALI) \
			$(PARTS_TRANSLATOR) \
			$(PARTS_SUBSTMATRIX) \
			$(PARTS_PHYLOLIB)

PARTS_PHYLOLIB= Tree.cpp HelpersTree.cpp \
                        ImplTree.cpp \
                        PhyloMatrix.cpp HelpersPhyloMatrix.cpp \
                        ImplPhyloMatrix.cpp ImplPhyloMatrixSymmetric.cpp \
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
			ImplNoWeightor.cpp ImplWeightor.cpp ImplWeightorHenikoff.cpp ImplWeightorHenikoffKimmen.cpp \
			ImplRegularizor.cpp  ImplNoRegularizor.cpp \
			ImplRegularizorDirichlet.cpp ImplRegularizorDirichletHash.cpp \
			ImplRegularizorDirichletInterpolate.cpp ImplRegularizorDirichletPrecomputed.cpp \
	                ImplLogOddor.cpp ImplLogOddorDirichlet.cpp ImplLogOddorUniform.cpp 

PARTS_ALIGNATA =      Alignata.cpp HelpersAlignata.cpp HelpersAlignataPrivate.cpp \
			ImplAlignata.cpp \
			ImplAlignataVector.cpp \
			ImplAlignataSet.cpp \
			ImplAlignataMatrix.cpp ImplAlignataMatrixRow.cpp \
			ImplAlignataMatrixDiagonal.cpp ImplAlignataMatrixUnsorted.cpp \
			ImplAlignataSetCol.cpp ImplAlignataHash.cpp ImplAlignataHashDiagonal.cpp 

PARTS_ALIGNATOR =	Alignator.cpp HelpersAlignator.cpp \
			ImplAlignator.cpp ImplAlignatorDP.cpp \
			ImplAlignatorDPFull.cpp \
			ImplAlignatorDots.cpp ImplAlignatorDotsSquared.cpp ImplAlignatorDotsSquaredDiagonal.cpp \
			ImplAlignatorIdentity.cpp ImplAlignatorSimilarity.cpp ImplAlignatorTuples.cpp \
			ImplAlignatorDummy.cpp \
			ImplAlignatorFragments.cpp ImplAlignatorFragmentsSquared.cpp \
			ImplAlignatorPublishAlignata.cpp

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
			ImplRenderer.cpp ImplRendererMView.cpp ImplRendererColumn.cpp \
			ImplAlignatum.cpp ImplAlignatumFasta.cpp 

PARTS_STATISTICS =    Statistics.cpp

PARTS_TRANSLATOR=	HelpersTranslator.cpp \
			ImplTranslator.cpp ImplTranslatorSequence.cpp \
			ImplTranslatorTCO.cpp ImplTranslatorSequenceUC.cpp

PARTS_SUBSTMATRIX=   	SubstitutionMatrix.cpp HelpersSubstitutionMatrix.cpp \
			ImplSubstitutionMatrix.cpp ImplSubstitutionMatrixAA.cpp

lib_LTLIBRARIES = libalignlib.la  

libalignlib_la_SOURCES = $(PARTS)

LIBOBJS= $(PARTS)


