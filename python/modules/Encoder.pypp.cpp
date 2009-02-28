// This file has been generated by Py++.

#include "boost/python.hpp"
#include "includes.h"
#include "iostream"
#include "cstdio"
#include "streambuf"
#include "Encoder.pypp.hpp"

namespace bp = boost::python;

class std_obuf: public std::streambuf 
        {
      public:
        std_obuf(std::FILE* file): m_file(file) {}
      protected:
        std::streambuf::int_type overflow(std::streambuf::int_type c) 
        {
          return std::fputc(c, m_file) ==EOF? std::streambuf::traits_type::eof(): c;
        }
      private:
        FILE* m_file;
      };

    class std_ibuf: public std::streambuf 
        {
      public:
          std_ibuf(std::FILE* file): m_file(file) {}
      protected:
          std::streambuf::int_type underflow() { 
           int c = std::getc(m_file);
           if (c != EOF) 
               std::ungetc(c, m_file);
            return c;
           }

          std::streambuf::int_type uflow() {
           return std::getc(m_file);
          }

          std::streambuf::int_type pbackfail(int c = EOF) {
            return c != EOF ? std::ungetc(c, m_file) : EOF;
          }
    private:
          FILE* m_file;      
          
      };

template<class T>
      void wrapper_for_save(const T & a, PyObject* fp) 
      {
        if (!PyFile_Check(fp)) 
        {
          throw boost::python::error_already_set();
        }
        std::FILE* f = PyFile_AsFile(fp);
        std_obuf buf(f);
        std::ostream os(&buf);
        a.save( os );
      }

void register_Encoder_class(){

    { //::alignlib::Encoder
        typedef bp::class_< alignlib::Encoder, bp::bases< alignlib::AlignlibBase >, boost::noncopyable > Encoder_exposer_t;
        Encoder_exposer_t Encoder_exposer = Encoder_exposer_t( "Encoder", bp::no_init );
        bp::scope Encoder_scope( Encoder_exposer );
        { //::alignlib::Encoder::countChars
        
            typedef ::alignlib::Position ( ::alignlib::Encoder::*countChars_function_type )( ::std::string const & ) const;
            
            Encoder_exposer.def( 
                "countChars"
                , countChars_function_type( &::alignlib::Encoder::countChars )
                , ( bp::arg("src") ) );
        
        }
        { //::alignlib::Encoder::decode
        
            typedef ::std::string ( ::alignlib::Encoder::*decode_function_type )( ::alignlib::ResidueVector const & ) const;
            
            Encoder_exposer.def( 
                "decode"
                , decode_function_type( &::alignlib::Encoder::decode )
                , ( bp::arg("src") ) );
        
        }
        { //::alignlib::Encoder::decode
        
            typedef char ( ::alignlib::Encoder::*decode_function_type )( ::alignlib::Residue const ) const;
            
            Encoder_exposer.def( 
                "decode"
                , decode_function_type( &::alignlib::Encoder::decode )
                , ( bp::arg("src") ) );
        
        }
        { //::alignlib::Encoder::encode
        
            typedef ::alignlib::ResidueVector ( ::alignlib::Encoder::*encode_function_type )( ::std::string const & ) const;
            
            Encoder_exposer.def( 
                "encode"
                , encode_function_type( &::alignlib::Encoder::encode )
                , ( bp::arg("src") ) );
        
        }
        { //::alignlib::Encoder::encode
        
            typedef ::alignlib::Residue ( ::alignlib::Encoder::*encode_function_type )( char const ) const;
            
            Encoder_exposer.def( 
                "encode"
                , encode_function_type( &::alignlib::Encoder::encode )
                , ( bp::arg("arg0") ) );
        
        }
        { //::alignlib::Encoder::getAlphabet
        
            typedef ::std::string ( ::alignlib::Encoder::*getAlphabet_function_type )(  ) const;
            
            Encoder_exposer.def( 
                "getAlphabet"
                , getAlphabet_function_type( &::alignlib::Encoder::getAlphabet ) );
        
        }
        { //::alignlib::Encoder::getAlphabetSize
        
            typedef int ( ::alignlib::Encoder::*getAlphabetSize_function_type )(  ) const;
            
            Encoder_exposer.def( 
                "getAlphabetSize"
                , getAlphabetSize_function_type( &::alignlib::Encoder::getAlphabetSize ) );
        
        }
        { //::alignlib::Encoder::getAlphabetType
        
            typedef ::alignlib::AlphabetType ( ::alignlib::Encoder::*getAlphabetType_function_type )(  ) const;
            
            Encoder_exposer.def( 
                "getAlphabetType"
                , getAlphabetType_function_type( &::alignlib::Encoder::getAlphabetType ) );
        
        }
        { //::alignlib::Encoder::getClone
        
            typedef ::alignlib::HEncoder ( ::alignlib::Encoder::*getClone_function_type )(  ) const;
            
            Encoder_exposer.def( 
                "getClone"
                , getClone_function_type( &::alignlib::Encoder::getClone ) );
        
        }
        { //::alignlib::Encoder::getGapChar
        
            typedef char ( ::alignlib::Encoder::*getGapChar_function_type )(  ) const;
            
            Encoder_exposer.def( 
                "getGapChar"
                , getGapChar_function_type( &::alignlib::Encoder::getGapChar ) );
        
        }
        { //::alignlib::Encoder::getGapChars
        
            typedef ::std::string ( ::alignlib::Encoder::*getGapChars_function_type )(  ) const;
            
            Encoder_exposer.def( 
                "getGapChars"
                , getGapChars_function_type( &::alignlib::Encoder::getGapChars ) );
        
        }
        { //::alignlib::Encoder::getGapCode
        
            typedef ::alignlib::Residue ( ::alignlib::Encoder::*getGapCode_function_type )(  ) const;
            
            Encoder_exposer.def( 
                "getGapCode"
                , getGapCode_function_type( &::alignlib::Encoder::getGapCode ) );
        
        }
        { //::alignlib::Encoder::getMap
        
            typedef ::alignlib::HResidueVector ( ::alignlib::Encoder::*getMap_function_type )( ::std::string const & ) const;
            
            Encoder_exposer.def( 
                "getMap"
                , getMap_function_type( &::alignlib::Encoder::getMap )
                , ( bp::arg("alphabet") ) );
        
        }
        { //::alignlib::Encoder::getMaskChar
        
            typedef char ( ::alignlib::Encoder::*getMaskChar_function_type )(  ) const;
            
            Encoder_exposer.def( 
                "getMaskChar"
                , getMaskChar_function_type( &::alignlib::Encoder::getMaskChar ) );
        
        }
        { //::alignlib::Encoder::getMaskChars
        
            typedef ::std::string ( ::alignlib::Encoder::*getMaskChars_function_type )(  ) const;
            
            Encoder_exposer.def( 
                "getMaskChars"
                , getMaskChars_function_type( &::alignlib::Encoder::getMaskChars ) );
        
        }
        { //::alignlib::Encoder::getMaskCode
        
            typedef ::alignlib::Residue ( ::alignlib::Encoder::*getMaskCode_function_type )(  ) const;
            
            Encoder_exposer.def( 
                "getMaskCode"
                , getMaskCode_function_type( &::alignlib::Encoder::getMaskCode ) );
        
        }
        { //::alignlib::Encoder::getNew
        
            typedef ::alignlib::HEncoder ( ::alignlib::Encoder::*getNew_function_type )(  ) const;
            
            Encoder_exposer.def( 
                "getNew"
                , getNew_function_type( &::alignlib::Encoder::getNew ) );
        
        }
        { //::alignlib::Encoder::isValidChar
        
            typedef bool ( ::alignlib::Encoder::*isValidChar_function_type )( char const ) const;
            
            Encoder_exposer.def( 
                "isValidChar"
                , isValidChar_function_type( &::alignlib::Encoder::isValidChar )
                , ( bp::arg("c") ) );
        
        }
        { //::alignlib::Encoder::map
        
            typedef ::alignlib::HResidueVector ( ::alignlib::Encoder::*map_function_type )( ::alignlib::HEncoder const & ) const;
            
            Encoder_exposer.def( 
                "map"
                , map_function_type( &::alignlib::Encoder::map )
                , ( bp::arg("other") ) );
        
        }
        { //::alignlib::Encoder::write
        
            typedef void ( ::alignlib::Encoder::*write_function_type )( ::std::ostream & ) const;
            
            Encoder_exposer.def( 
                "write"
                , write_function_type( &::alignlib::Encoder::write )
                , ( bp::arg("arg0") ) );
        
        }
        Encoder_exposer.def( bp::self_ns::str( bp::self ) );
        bp::register_ptr_to_python< boost::shared_ptr< alignlib::Encoder > >();
        bp::implicitly_convertible< boost::shared_ptr< alignlib::Encoder >, boost::shared_ptr< alignlib::AlignlibBase > >();
        Encoder_exposer.def( "save", wrapper_for_save<alignlib::Encoder> );
    }

}
