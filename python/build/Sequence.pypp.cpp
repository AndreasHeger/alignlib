// This file has been generated by Py++.

#include "boost/python.hpp"
#include "includes.h"
#include "iostream"
#include "cstdio"
#include "Sequence.pypp.hpp"

namespace bp = boost::python;

void register_Sequence_class(){

    { //::alignlib::Sequence
        typedef bp::class_< alignlib::Sequence, bp::bases< alignlib::Alignandum >, boost::noncopyable > Sequence_exposer_t;
        Sequence_exposer_t Sequence_exposer = Sequence_exposer_t( "Sequence", bp::no_init );
        bp::scope Sequence_scope( Sequence_exposer );
        Sequence_exposer.def( bp::self_ns::str( bp::self ) );
        bp::register_ptr_to_python< boost::shared_ptr< alignlib::Sequence > >();
        bp::implicitly_convertible< boost::shared_ptr< alignlib::Sequence >, boost::shared_ptr< alignlib::Alignandum > >();
    }

}