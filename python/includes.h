#include "alignlib.h"
#include "Alignandum.h"
#include "Alignata.h"
#include "AlignataIterator.h"
#include "Alignator.h"
#include "Alignatum.h"
#include "AlignException.h"
#include "AlignlibDebug.h"
#include "AlignlibMethods.h"
// #include "Dottor.h"
#include "Fragmentor.h"
#include "HelpersAlignandum.h"
#include "HelpersAlignata.h"
#include "HelpersAlignator.h"
#include "HelpersAlignatum.h"
#include "HelpersFragmentor.h"
#include "HelpersIterator2D.h"
#include "HelpersLogOddor.h"
#include "HelpersMultipleAlignment.h"
#include "HelpersProfile.h"
#include "HelpersRegularizor.h"
#include "HelpersRenderer.h"
#include "HelpersScorer.h"
#include "HelpersSequence.h"
#include "HelpersSubstitutionMatrix.h"
#include "HelpersTranslator.h"
#include "HelpersWeightor.h"
#include "Iterator2DFull.h"
#include "Iterator2D.h"
#include "LogOddor.h"
#include "Matrix.h"
#include "MultipleAlignment.h"
#include "Regularizor.h"
#include "Renderer.h"
#include "Scorer.h"
#include "Statistics.h"
#include "SubstitutionMatrix.h"
#include "Translator.h"
#include "Weightor.h"

namespace py_details
{
  inline void instantiate() 
  {
    alignlib::Matrix<int> m1(1,1);
    alignlib::Matrix<double> m2(1,1);
    alignlib::Matrix<unsigned int> m3(1,1);
  }
}
