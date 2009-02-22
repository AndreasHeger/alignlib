// This file has been generated by Py++.

#include "boost/python.hpp"
#include "includes.h"
#include "iostream"
#include "cstdio"
#include "Iterator2D.pypp.hpp"

namespace bp = boost::python;

void register_Iterator2D_class(){

    { //::alignlib::Iterator2D
        typedef bp::class_< alignlib::Iterator2D, boost::noncopyable > Iterator2D_exposer_t;
        Iterator2D_exposer_t Iterator2D_exposer = Iterator2D_exposer_t( "Iterator2D", bp::no_init );
        bp::scope Iterator2D_scope( Iterator2D_exposer );
        { //::alignlib::Iterator2D::col_back
        
            typedef ::alignlib::Position ( ::alignlib::Iterator2D::*col_back_function_type )( ::alignlib::Position ) const;
            
            Iterator2D_exposer.def( 
                "col_back"
                , col_back_function_type( &::alignlib::Iterator2D::col_back )
                , ( bp::arg("row")=(int)(-0x00000000000000001) ) );
        
        }
        { //::alignlib::Iterator2D::col_begin
        
            typedef ::alignlib::const_countable_iterator< int > ( ::alignlib::Iterator2D::*col_begin_function_type )( ::alignlib::Position ) const;
            
            Iterator2D_exposer.def( 
                "col_begin"
                , col_begin_function_type( &::alignlib::Iterator2D::col_begin )
                , ( bp::arg("row")=(int)(-0x00000000000000001) ) );
        
        }
        { //::alignlib::Iterator2D::col_end
        
            typedef ::alignlib::const_countable_iterator< int > ( ::alignlib::Iterator2D::*col_end_function_type )( ::alignlib::Position ) const;
            
            Iterator2D_exposer.def( 
                "col_end"
                , col_end_function_type( &::alignlib::Iterator2D::col_end )
                , ( bp::arg("row")=(int)(-0x00000000000000001) ) );
        
        }
        { //::alignlib::Iterator2D::col_front
        
            typedef ::alignlib::Position ( ::alignlib::Iterator2D::*col_front_function_type )( ::alignlib::Position ) const;
            
            Iterator2D_exposer.def( 
                "col_front"
                , col_front_function_type( &::alignlib::Iterator2D::col_front )
                , ( bp::arg("row")=(int)(-0x00000000000000001) ) );
        
        }
        { //::alignlib::Iterator2D::col_size
        
            typedef ::alignlib::Position ( ::alignlib::Iterator2D::*col_size_function_type )( ::alignlib::Position ) const;
            
            Iterator2D_exposer.def( 
                "col_size"
                , col_size_function_type( &::alignlib::Iterator2D::col_size )
                , ( bp::arg("row")=(int)(-0x00000000000000001) ) );
        
        }
        { //::alignlib::Iterator2D::getClone
        
            typedef ::alignlib::HIterator2D ( ::alignlib::Iterator2D::*getClone_function_type )(  ) const;
            
            Iterator2D_exposer.def( 
                "getClone"
                , getClone_function_type( &::alignlib::Iterator2D::getClone ) );
        
        }
        { //::alignlib::Iterator2D::getNew
        
            typedef ::alignlib::HIterator2D ( ::alignlib::Iterator2D::*getNew_function_type )(  ) const;
            
            Iterator2D_exposer.def( 
                "getNew"
                , getNew_function_type( &::alignlib::Iterator2D::getNew ) );
        
        }
        { //::alignlib::Iterator2D::getNew
        
            typedef ::alignlib::HIterator2D ( ::alignlib::Iterator2D::*getNew_function_type )( ::alignlib::HAlignandum const &,::alignlib::HAlignandum const & ) const;
            
            Iterator2D_exposer.def( 
                "getNew"
                , getNew_function_type( &::alignlib::Iterator2D::getNew )
                , ( bp::arg("row"), bp::arg("col") ) );
        
        }
        { //::alignlib::Iterator2D::resetRanges
        
            typedef void ( ::alignlib::Iterator2D::*resetRanges_function_type )( ::alignlib::HAlignandum const &,::alignlib::HAlignandum const & ) ;
            
            Iterator2D_exposer.def( 
                "resetRanges"
                , resetRanges_function_type( &::alignlib::Iterator2D::resetRanges )
                , ( bp::arg("row"), bp::arg("col") ) );
        
        }
        { //::alignlib::Iterator2D::row_back
        
            typedef ::alignlib::Position ( ::alignlib::Iterator2D::*row_back_function_type )( ::alignlib::Position ) const;
            
            Iterator2D_exposer.def( 
                "row_back"
                , row_back_function_type( &::alignlib::Iterator2D::row_back )
                , ( bp::arg("col")=(int)(-0x00000000000000001) ) );
        
        }
        { //::alignlib::Iterator2D::row_begin
        
            typedef ::alignlib::const_countable_iterator< int > ( ::alignlib::Iterator2D::*row_begin_function_type )( ::alignlib::Position ) const;
            
            Iterator2D_exposer.def( 
                "row_begin"
                , row_begin_function_type( &::alignlib::Iterator2D::row_begin )
                , ( bp::arg("col")=(int)(-0x00000000000000001) ) );
        
        }
        { //::alignlib::Iterator2D::row_end
        
            typedef ::alignlib::const_countable_iterator< int > ( ::alignlib::Iterator2D::*row_end_function_type )( ::alignlib::Position ) const;
            
            Iterator2D_exposer.def( 
                "row_end"
                , row_end_function_type( &::alignlib::Iterator2D::row_end )
                , ( bp::arg("col")=(int)(-0x00000000000000001) ) );
        
        }
        { //::alignlib::Iterator2D::row_front
        
            typedef ::alignlib::Position ( ::alignlib::Iterator2D::*row_front_function_type )( ::alignlib::Position ) const;
            
            Iterator2D_exposer.def( 
                "row_front"
                , row_front_function_type( &::alignlib::Iterator2D::row_front )
                , ( bp::arg("col")=(int)(-0x00000000000000001) ) );
        
        }
        { //::alignlib::Iterator2D::row_size
        
            typedef ::alignlib::Position ( ::alignlib::Iterator2D::*row_size_function_type )( ::alignlib::Position ) const;
            
            Iterator2D_exposer.def( 
                "row_size"
                , row_size_function_type( &::alignlib::Iterator2D::row_size )
                , ( bp::arg("col")=(int)(-0x00000000000000001) ) );
        
        }
        bp::register_ptr_to_python< boost::shared_ptr< alignlib::Iterator2D > >();
    }

}
