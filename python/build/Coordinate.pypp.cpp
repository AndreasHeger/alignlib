// This file has been generated by Py++.

#include "boost/python.hpp"
#include "includes.h"
#include "iostream"
#include "cstdio"
#include "Coordinate.pypp.hpp"

namespace bp = boost::python;

void register_Coordinate_class(){

    { //::alignlib::Coordinate
        typedef bp::class_< alignlib::Coordinate > Coordinate_exposer_t;
        Coordinate_exposer_t Coordinate_exposer = Coordinate_exposer_t( "Coordinate" );
        bp::scope Coordinate_scope( Coordinate_exposer );
        Coordinate_exposer.def_readwrite( "col", &alignlib::Coordinate::col );
        Coordinate_exposer.def_readwrite( "row", &alignlib::Coordinate::row );
    }

}