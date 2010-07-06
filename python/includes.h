
#include "alignlib.h"

#include <iostream>

namespace py_details
{
  inline void instantiate()
  {
    // call these typedefs to register them for alias creation
    // otherwise: ugly alias warnings will appear
    size_t t;
    t = sizeof(alignlib::Matrix<int>);
    t = sizeof(alignlib::Matrix<double>);
    t = sizeof(alignlib::Matrix<unsigned int>);
    /* not needed any more? */
    /* t = sizeof(std::vector<int>); */
    /* t = sizeof(std::vector<unsigned int>); */
    /* t = sizeof(std::vector<long>); */
    /* t = sizeof(std::vector<unsigned long>); */
    /* t = sizeof(alignlib::SubstitutionMatrix); */
    /* t = sizeof(alignlib::MutationMatrix); */
    /* t = sizeof(alignlib::FragmentVector); */
    /* t = sizeof(alignlib::PositionVector); */
    /* t = sizeof(alignlib::ScoreVector); */
    /* t = sizeof(alignlib::EntropyVector); */
    /* t = sizeof(alignlib::CountVector); */
    /* t = sizeof(alignlib::SequenceWeights); */
    /* t = sizeof(alignlib::NodeVector); */
    /* t = sizeof(alignlib::StringVector); */
    /* t = sizeof(alignlib::HAlignandum); */
    /*   t = sizeof(alignlib::HAlignment); */
    /*   t = sizeof(alignlib::HAlignator); */
    /*   t = sizeof(alignlib::HMultipleAlignment); */
    /*   t = sizeof(boost::shared_ptr<alignlib::Alignandum>); */
    /*   t = sizeof(boost::shared_ptr<alignlib::Alignment>); */
    /*   t = sizeof(boost::shared_ptr<alignlib::MultipleAlignment>); */
    /*   t = sizeof(boost::shared_ptr<alignlib::Alignator>); */
    /*   t = sizeof(boost::shared_ptr<alignlib::Alignment>); */


  }
}


