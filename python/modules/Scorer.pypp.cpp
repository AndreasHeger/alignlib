// This file has been generated by Py++.

#include "boost/python.hpp"
#include "includes.h"
#include "iostream"
#include "cstdio"
#include "Scorer.pypp.hpp"

namespace bp = boost::python;

void register_Scorer_class(){

    { //::alignlib::Scorer
        typedef bp::class_< alignlib::Scorer, bp::bases< alignlib::AlignlibBase >, boost::noncopyable > Scorer_exposer_t;
        Scorer_exposer_t Scorer_exposer = Scorer_exposer_t( "Scorer", bp::no_init );
        bp::scope Scorer_scope( Scorer_exposer );
        { //::alignlib::Scorer::getClone
        
            typedef ::alignlib::HScorer ( ::alignlib::Scorer::*getClone_function_type )(  ) const;
            
            Scorer_exposer.def( 
                "getClone"
                , getClone_function_type( &::alignlib::Scorer::getClone ) );
        
        }
        { //::alignlib::Scorer::getNew
        
            typedef ::alignlib::HScorer ( ::alignlib::Scorer::*getNew_function_type )(  ) const;
            
            Scorer_exposer.def( 
                "getNew"
                , getNew_function_type( &::alignlib::Scorer::getNew ) );
        
        }
        { //::alignlib::Scorer::getNew
        
            typedef ::alignlib::HScorer ( ::alignlib::Scorer::*getNew_function_type )( ::alignlib::HAlignandum const &,::alignlib::HAlignandum const & ) const;
            
            Scorer_exposer.def( 
                "getNew"
                , getNew_function_type( &::alignlib::Scorer::getNew )
                , ( bp::arg("row"), bp::arg("col") ) );
        
        }
        { //::alignlib::Scorer::getScore
        
            typedef ::alignlib::Score ( ::alignlib::Scorer::*getScore_function_type )( ::alignlib::Position const &,::alignlib::Position const & ) const;
            
            Scorer_exposer.def( 
                "getScore"
                , getScore_function_type( &::alignlib::Scorer::getScore )
                , ( bp::arg("row"), bp::arg("col") ) );
        
        }
        bp::register_ptr_to_python< boost::shared_ptr< alignlib::Scorer > >();
        bp::implicitly_convertible< boost::shared_ptr< alignlib::Scorer >, boost::shared_ptr< alignlib::AlignlibBase > >();
    }

}
