#include "alignlib.h"

#include "Alignandum.h"
#include "HelpersAlignandum.h"
#include "HelpersProfile.h"
#include "HelpersSequence.h"

#include "Alignment.h"
#include "HelpersAlignment.h"

#include "AlignmentIterator.h"

#include "Alignator.h"
#include "HelpersAlignator.h"

#include "Alignatum.h"
#include "HelpersAlignatum.h"

#include "AlignException.h"
#include "AlignlibDebug.h"
#include "AlignlibMethods.h"

#include "Fragmentor.h"
#include "HelpersFragmentor.h"

#include "Iterator2D.h"
#include "Iterator2DFull.h"
#include "HelpersIterator2D.h"

#include "LogOddor.h"
#include "HelpersLogOddor.h"

#include "Matrix.h"

#include "MultipleAlignment.h"
#include "HelpersMultipleAlignment.h"

#include "Regularizor.h"
#include "HelpersRegularizor.h"

#include "Renderer.h"
#include "HelpersRenderer.h"

#include "Scorer.h"
#include "HelpersScorer.h"

#include "Statistics.h"

#include "SubstitutionMatrix.h"
#include "HelpersSubstitutionMatrix.h"

#include "Translator.h"
#include "HelpersTranslator.h"

#include "Weightor.h"
#include "HelpersWeightor.h"

#include "Distor.h"
#include "HelpersDistor.h"

#include "Treetor.h"
#include "HelpersTreetor.h"

#include "PhyloMatrix.h"
#include "HelpersPhyloMatrix.h"

#include "Tree.h"
#include "HelpersTree.h"

#include "AlignedBlocks.h"

#include <iostream>

namespace py_details
{
  inline void instantiate() 
  {
    sizeof(alignlib::Matrix<int>);
    sizeof(alignlib::Matrix<double>);
    sizeof(alignlib::Matrix<unsigned int>);
  }
}