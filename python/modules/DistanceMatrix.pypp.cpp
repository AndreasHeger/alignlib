// This file has been generated by Py++.

#include "boost/python.hpp"
#include "includes.h"
#include "iostream"
#include "cstdio"
#include "DistanceMatrix.pypp.hpp"

namespace bp = boost::python;

void register_DistanceMatrix_class(){

    { //::alignlib::DistanceMatrix
        typedef bp::class_< alignlib::DistanceMatrix, bp::bases< alignlib::AlignlibBase >, boost::noncopyable > DistanceMatrix_exposer_t;
        DistanceMatrix_exposer_t DistanceMatrix_exposer = DistanceMatrix_exposer_t( "DistanceMatrix", bp::no_init );
        bp::scope DistanceMatrix_scope( DistanceMatrix_exposer );
        { //::alignlib::DistanceMatrix::getClone
        
            typedef ::alignlib::HDistanceMatrix ( ::alignlib::DistanceMatrix::*getClone_function_type )(  ) const;
            
            DistanceMatrix_exposer.def( 
                "getClone"
                , getClone_function_type( &::alignlib::DistanceMatrix::getClone ) );
        
        }
        { //::alignlib::DistanceMatrix::getElement
        
            typedef ::alignlib::DistanceMatrixValue ( ::alignlib::DistanceMatrix::*getElement_function_type )( ::alignlib::DistanceMatrixSize,::alignlib::DistanceMatrixSize ) const;
            
            DistanceMatrix_exposer.def( 
                "getElement"
                , getElement_function_type( &::alignlib::DistanceMatrix::getElement )
                , ( bp::arg("row"), bp::arg("col") ) );
        
        }
        { //::alignlib::DistanceMatrix::getMaximum
        
            typedef ::alignlib::DistanceMatrixValue ( ::alignlib::DistanceMatrix::*getMaximum_function_type )(  ) const;
            
            DistanceMatrix_exposer.def( 
                "getMaximum"
                , getMaximum_function_type( &::alignlib::DistanceMatrix::getMaximum ) );
        
        }
        { //::alignlib::DistanceMatrix::getMaximum
        
            typedef ::alignlib::DistanceMatrixValue ( ::alignlib::DistanceMatrix::*getMaximum_function_type )( ::alignlib::Coordinate & ) const;
            
            DistanceMatrix_exposer.def( 
                "getMaximum"
                , getMaximum_function_type( &::alignlib::DistanceMatrix::getMaximum )
                , ( bp::arg("dest") ) );
        
        }
        { //::alignlib::DistanceMatrix::getMinimum
        
            typedef ::alignlib::DistanceMatrixValue ( ::alignlib::DistanceMatrix::*getMinimum_function_type )(  ) const;
            
            DistanceMatrix_exposer.def( 
                "getMinimum"
                , getMinimum_function_type( &::alignlib::DistanceMatrix::getMinimum ) );
        
        }
        { //::alignlib::DistanceMatrix::getMinimum
        
            typedef ::alignlib::DistanceMatrixValue ( ::alignlib::DistanceMatrix::*getMinimum_function_type )( ::alignlib::Coordinate & ) const;
            
            DistanceMatrix_exposer.def( 
                "getMinimum"
                , getMinimum_function_type( &::alignlib::DistanceMatrix::getMinimum )
                , ( bp::arg("dest") ) );
        
        }
        { //::alignlib::DistanceMatrix::getNew
        
            typedef ::alignlib::HDistanceMatrix ( ::alignlib::DistanceMatrix::*getNew_function_type )(  ) const;
            
            DistanceMatrix_exposer.def( 
                "getNew"
                , getNew_function_type( &::alignlib::DistanceMatrix::getNew ) );
        
        }
        { //::alignlib::DistanceMatrix::getSize
        
            typedef ::alignlib::DistanceMatrixSize ( ::alignlib::DistanceMatrix::*getSize_function_type )(  ) const;
            
            DistanceMatrix_exposer.def( 
                "getSize"
                , getSize_function_type( &::alignlib::DistanceMatrix::getSize ) );
        
        }
        { //::alignlib::DistanceMatrix::getWidth
        
            typedef ::alignlib::DistanceMatrixSize ( ::alignlib::DistanceMatrix::*getWidth_function_type )(  ) const;
            
            DistanceMatrix_exposer.def( 
                "getWidth"
                , getWidth_function_type( &::alignlib::DistanceMatrix::getWidth ) );
        
        }
        { //::alignlib::DistanceMatrix::read
        
            typedef void ( ::alignlib::DistanceMatrix::*read_function_type )( ::std::istream & ) const;
            
            DistanceMatrix_exposer.def( 
                "read"
                , read_function_type( &::alignlib::DistanceMatrix::read )
                , ( bp::arg("input") ) );
        
        }
        { //::alignlib::DistanceMatrix::setElement
        
            typedef void ( ::alignlib::DistanceMatrix::*setElement_function_type )( ::alignlib::DistanceMatrixSize,::alignlib::DistanceMatrixSize,::alignlib::DistanceMatrixValue ) ;
            
            DistanceMatrix_exposer.def( 
                "setElement"
                , setElement_function_type( &::alignlib::DistanceMatrix::setElement )
                , ( bp::arg("row"), bp::arg("col"), bp::arg("value") ) );
        
        }
        { //::alignlib::DistanceMatrix::setWidth
        
            typedef void ( ::alignlib::DistanceMatrix::*setWidth_function_type )( ::alignlib::DistanceMatrixSize ) ;
            
            DistanceMatrix_exposer.def( 
                "setWidth"
                , setWidth_function_type( &::alignlib::DistanceMatrix::setWidth )
                , ( bp::arg("width") ) );
        
        }
        { //::alignlib::DistanceMatrix::shrink
        
            typedef void ( ::alignlib::DistanceMatrix::*shrink_function_type )(  ) ;
            
            DistanceMatrix_exposer.def( 
                "shrink"
                , shrink_function_type( &::alignlib::DistanceMatrix::shrink ) );
        
        }
        { //::alignlib::DistanceMatrix::swap
        
            typedef void ( ::alignlib::DistanceMatrix::*swap_function_type )( ::alignlib::DistanceMatrixSize,::alignlib::DistanceMatrixSize ) ;
            
            DistanceMatrix_exposer.def( 
                "swap"
                , swap_function_type( &::alignlib::DistanceMatrix::swap )
                , ( bp::arg("a"), bp::arg("b") ) );
        
        }
        { //::alignlib::DistanceMatrix::write
        
            typedef void ( ::alignlib::DistanceMatrix::*write_function_type )( ::std::ostream & ) const;
            
            DistanceMatrix_exposer.def( 
                "write"
                , write_function_type( &::alignlib::DistanceMatrix::write )
                , ( bp::arg("output") ) );
        
        }
        DistanceMatrix_exposer.def( bp::self_ns::str( bp::self ) );
        bp::register_ptr_to_python< boost::shared_ptr< alignlib::DistanceMatrix > >();
        bp::implicitly_convertible< boost::shared_ptr< alignlib::DistanceMatrix >, boost::shared_ptr< alignlib::AlignlibBase > >();
    }

}
