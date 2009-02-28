// This file has been generated by Py++.

#include "boost/python.hpp"
#include "includes.h"
#include "iostream"
#include "cstdio"
#include "Distor.pypp.hpp"

namespace bp = boost::python;

void register_Distor_class(){

    { //::alignlib::Distor
        typedef bp::class_< alignlib::Distor, bp::bases< alignlib::AlignlibBase >, boost::noncopyable > Distor_exposer_t;
        Distor_exposer_t Distor_exposer = Distor_exposer_t( "Distor", bp::no_init );
        bp::scope Distor_scope( Distor_exposer );
        { //::alignlib::Distor::calculateDistance
        
            typedef ::alignlib::DistanceMatrixValue ( ::alignlib::Distor::*calculateDistance_function_type )( ::std::string const &,::std::string const & ) const;
            
            Distor_exposer.def( 
                "calculateDistance"
                , calculateDistance_function_type( &::alignlib::Distor::calculateDistance )
                , ( bp::arg("s_row_1"), bp::arg("s_row_2") ) );
        
        }
        { //::alignlib::Distor::calculateMatrix
        
            typedef void ( ::alignlib::Distor::*calculateMatrix_function_type )( ::alignlib::HDistanceMatrix &,::alignlib::HMultipleAlignment const & ) const;
            
            Distor_exposer.def( 
                "calculateMatrix"
                , calculateMatrix_function_type( &::alignlib::Distor::calculateMatrix )
                , ( bp::arg("dest"), bp::arg("mali") ) );
        
        }
        { //::alignlib::Distor::getClone
        
            typedef ::alignlib::HDistor ( ::alignlib::Distor::*getClone_function_type )(  ) const;
            
            Distor_exposer.def( 
                "getClone"
                , getClone_function_type( &::alignlib::Distor::getClone ) );
        
        }
        { //::alignlib::Distor::getMaximumPossibleDistance
        
            typedef ::alignlib::DistanceMatrixValue ( ::alignlib::Distor::*getMaximumPossibleDistance_function_type )(  ) const;
            
            Distor_exposer.def( 
                "getMaximumPossibleDistance"
                , getMaximumPossibleDistance_function_type( &::alignlib::Distor::getMaximumPossibleDistance ) );
        
        }
        { //::alignlib::Distor::getNew
        
            typedef ::alignlib::HDistor ( ::alignlib::Distor::*getNew_function_type )(  ) const;
            
            Distor_exposer.def( 
                "getNew"
                , getNew_function_type( &::alignlib::Distor::getNew ) );
        
        }
        Distor_exposer.def( bp::self_ns::str( bp::self ) );
        bp::register_ptr_to_python< boost::shared_ptr< alignlib::Distor > >();
        bp::implicitly_convertible< boost::shared_ptr< alignlib::Distor >, boost::shared_ptr< alignlib::AlignlibBase > >();
    }

}
