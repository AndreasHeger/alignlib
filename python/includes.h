#include "alignlib/alignlib.h"
#include "alignlib/Alignandum.h"
#include "alignlib/Alignata.h"
#include "alignlib/AlignataIterator.h"
#include "alignlib/Alignator.h"
#include "alignlib/Alignatum.h"
#include "alignlib/AlignException.h"
#include "alignlib/AlignlibDebug.h"
#include "alignlib/alignlib.h"
#include "alignlib/AlignlibMethods.h"
// #include "alignlib/Dottor.h"
#include "alignlib/Fragmentor.h"
#include "alignlib/HelpersAlignandum.h"
#include "alignlib/HelpersAlignata.h"
#include "alignlib/HelpersAlignator.h"
#include "alignlib/HelpersAlignatum.h"
#include "alignlib/HelpersFragmentor.h"
#include "alignlib/HelpersIterator2D.h"
#include "alignlib/HelpersLogOddor.h"
#include "alignlib/HelpersMultipleAlignment.h"
#include "alignlib/HelpersProfile.h"
#include "alignlib/HelpersRegularizor.h"
#include "alignlib/HelpersRenderer.h"
#include "alignlib/HelpersScorer.h"
#include "alignlib/HelpersSequence.h"
#include "alignlib/HelpersSubstitutionMatrix.h"
#include "alignlib/HelpersTranslator.h"
#include "alignlib/HelpersWeightor.h"
#include "alignlib/Iterator2DFull.h"
#include "alignlib/Iterator2D.h"
#include "alignlib/LogOddor.h"
#include "alignlib/Matrix.h"
#include "alignlib/MultipleAlignment.h"
#include "alignlib/Regularizor.h"
#include "alignlib/Renderer.h"
#include "alignlib/Scorer.h"
#include "alignlib/Statistics.h"
#include "alignlib/SubstitutionMatrix.h"
#include "alignlib/Translator.h"
#include "alignlib/Weightor.h"

namespace py_details
{
  inline void instantiate() 
  {
    alignlib::Matrix<int> m1(1,1);
    alignlib::Matrix<double> m2(1,1);
    alignlib::Matrix<unsigned int> m3(1,1);
  }
}
