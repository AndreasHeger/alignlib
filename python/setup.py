
USAGE = """python setup.py [options] command [command [...]]

The following commands are available:

* build: build and compile python extension to alignlib
* test: run some tests
* install: install the extension
"""

import re, sys, os, optparse, subprocess

from pyplusplus import module_builder, messages, function_transformers
from pyplusplus.module_builder import call_policies
from pyplusplus.decl_wrappers import \
     return_value_policy, manage_new_object, copy_const_reference, reference_existing_object, \
     return_self, return_arg

from pygccxml import declarations

def findBoost( options ):
    """find location of boost."""
    if not options.boost_dir:
        if "BOOST_ROOT" in os.environ:
            options.boost_dir = os.environ['BOOST_ROOT']
            return 

        for x in ("/usr/local/boost", "/opt/boost" ):
            if os.path.exists( x ):
               options.boost_dir = x 
               return

        raise "could not find BOOST. Please specify location of BOOST as option or set BOOST_ROOT environment variable."

def checkRequisites( options ):
    """check if boost is present."""
    
    nerrors = 0
    if not os.path.exists(options.boost_dir):
        nerrors += 1
        print "could not find boost directory %s" % options.boost_dir
                    
    if not os.path.exists(options.alignlib_lib_dir ):                
        nerrors += 1
        print "could not find library directory %s" % options.alignlib_lib_dir
    elif not os.path.exists( options.alignlib_lib_dir + "/libalignlib.so" ):
        nerrors += 1
        print "could not find alinglib shared library %s/libalignlib.so" % options.alignlib_lib_dir
        
    return nerrors

def addStreamBufClasses( mb ):
    """add wrappers to export mapping from c++ streams to python file objects."""
    declaration_code = \
    """
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
            
    """
    mb.add_declaration_code ( declaration_code, tail=True )

def exportSave( mb, options):
    """export save method in classes Alignandum."""
    
    declaration_code = \
"""
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
    """
    
    cls = mb.class_("Alignandum")
    cls.include_files.append( "streambuf" )
    registration_code = 'def( "save", wrapper_for_save<alignlib::%s> )' % cls.name
    cls.member_function( "save" ).exclude()
    mb.add_declaration_code( declaration_code, tail=True )
    cls.add_registration_code( registration_code )

def exportLoad( mb, options ):
    
    declaration_code = \
    """        
      alignlib::Alignandum * wrapper_for_load( PyObject * fp )
      {
          if (!PyFile_Check(fp))
          {
              throw boost::python::error_already_set();
          }
         std::FILE * f = PyFile_AsFile(fp);   
         std_ibuf buf(f);
         std::istream is(&buf);
         alignlib::Alignandum * a = alignlib::loadAlignandum( is );
         return a;
      }       
"""
    registration_code = 'def( "bp::loadAlignandum", wrapper_for_load, bp::return_value_policy< bp::manage_new_object >() );'     
    
    fun = mb.free_function( "loadAlignandum")
    fun.include_files.append( "streambuf" )
    fun.exclude()
    mb.add_declaration_code( declaration_code, tail = True )
    mb.add_registration_code( registration_code )

def export_writePairAlignment( mb ):
    """export writePairAlignment.
    
    This function has to map ostream to a python file.
    """
    declaration_code = \
"""
      void wrapper_for_writePairAlignment( PyObject* fp, 
          const alignlib::Alignandum * row, 
          const alignlib::Alignandum * col,
          const alignlib::Alignment * alignment ) 
      {
        if (!PyFile_Check(fp)) 
        {
          throw boost::python::error_already_set();
        }
        std::FILE* f = PyFile_AsFile(fp);
        std_obuf buf(f);
        std::ostream os(&buf);
        writePairAlignment( os, row, col, alignment );
      }
    """
    
    fun = mb.free_function("writePairAlignment")
    fun.include_files.append( "streambuf" )
    fun.exclude()
    registration_code = """bp::def( "writePairAlignment", 
        wrapper_for_writePairAlignment, 
        (bp::arg("output"), bp::arg("row"), bp::arg("col"), bp::arg("ali") ));"""
         
    mb.add_declaration_code( declaration_code, tail = True )
    mb.add_registration_code( registration_code )

def export_writeAlignmentCompressed( mb ):
    """export writePairAlignment.
    
    This function has to map ostream to a python file.
    """
    declaration_code = \
"""
      void wrapper_for_writeAlignmentCompressed( PyObject* fp, 
          const alignlib::Alignment * alignment )
      {
        if (!PyFile_Check(fp)) 
        {
          throw boost::python::error_already_set();
        }
        std::FILE* f = PyFile_AsFile(fp);
        std_obuf buf(f);
        std::ostream os(&buf);
        writeAlignmentCompressed( os, alignment );
      }
    """
    
    fun = mb.free_function("writeAlignmentCompressed")
    fun.include_files.append( "streambuf" )
    fun.exclude()
    registration_code = """bp::def( "writeAlignmentCompressed", 
        wrapper_for_writeAlignmentCompressed, 
        (bp::arg("output"), bp::arg("ali") ));"""
         
    mb.add_declaration_code( declaration_code, tail = True )
    mb.add_registration_code( registration_code )


def export_functions( mb ):
    """export utility functions."""
    ## include all membership functions
    mb.namespace("alignlib").free_functions().include()

    # set call policies for functions that return the same object
    # from the first argument
    for prefix in ("fill",
                   "add",
                   "substitute",
                   "reset",
                   "copy",
                   "combine",
                   "complement",
                   "rescore",
                   "flatten",
                   "filter",
                   "readAlignmentPairs",
                   "calculateAffineScore",
#                       "extractMultipleAlignment",                       
                   "rescaleProfileCounts",
                   "normalizeProfileCounts",
                   ):
        try:
            mb.free_functions( lambda mem_fun: mem_fun.name.startswith( prefix )).call_policies = return_self()
        except RuntimeError:
            sys.stderr.write("# could not find free function starting with %s\n" % prefix )

    # set call policies for functions that return get a default object
    # in this case the caller is not a new object.
    for prefix in ("getDefault", ):
        mb.free_functions( lambda mem_fun: mem_fun.name.startswith( prefix )).call_policies = \
                           return_value_policy( reference_existing_object )

    # Set call policies for functions that set a default object
    # The caller takes ownership of the old default object and
    # and passes control of the new default object to the library
    # This only worked for non-const argument types and not returning
    # an new object.
    for prefix in ("setDefault", ):
        
        for fun in mb.free_functions( lambda mem_fun: mem_fun.name.startswith( prefix ) ):
            
            cpointee = fun.name[len(prefix):]
            cls_pointee = mb.class_( cpointee )
            cls_pointee.held_type = 'std::auto_ptr< %s >' % cls_pointee.decl_string
                    
            fname = fun.name
            # set alias to original function name, otherwise ugly names will be created
            fun.add_transformation( function_transformers.transfer_ownership( 0 ), 
                                    alias = fname )

            # the following did not work, thus changed setDefault to return void
            # fun.call_policies = return_value_policy( manage_new_object )

    # other functions that return new objects, including the factory functions
    for prefix in ("extract", "read", "make", "load", "exportProfileFrequencies", "splitAlignment" ): 
        try:
            mb.free_functions( lambda mem_fun: mem_fun.name.startswith( prefix )).call_policies = \
                               return_value_policy( manage_new_object )
        except RuntimeError:
            print "could not find any function with prefix %s" % prefix
            
    #######################################################################################
    #######################################################################################
    #######################################################################################
    ## patches to exclude problematic functions
    for prefix in ("getMapResidue", "makeRendererColumn", "makeAlignatorDotsWrap", "getDefaultPalette" ):
        try:
            mb.free_functions( lambda mem_fun: mem_fun.name.startswith( prefix )).exclude()
        except RuntimeError:
            print "could not find declaration for %s" % prefix
            
    ## can't export makeProfile(const std::string &, int)
    ## couldn't figure out what to mach const std::string & with
    ## tried: const std::string &, std::string const &, and more
    mb.free_functions( name='makeProfile', arg_types=[None, "int"] ).exclude()
    
    ## deal with functions that work on streams
    export_writePairAlignment( mb )
    export_writeAlignmentCompressed( mb )

   
def export_classes( mb ):
    """export classes.
    
    These classes can be instantiated directly from python.
    """
    classes_to_export = set( ['AlignedBlocks',] )
    
    ## include all classes
    mb.classes( lambda x: x.name in classes_to_export ).include()

    mb.member_functions( lambda x: x.name == "copy" ).call_policies = \
        return_value_policy( reference_existing_object )

        
def export_interface_classes( mb ):
    """export virtual classes.
    
    These classes can not instantiated directly from python.
    """
    classes_to_export = set( ['Alignandum',
                              'MultipleAlignment',
                              'Alignator',
                              'Iterator',
    #                          'Dottor',
    #                          'Translator',
                              'SubstitutionMatrix',
                              'Fragmentor',
                              'Alignment',
                              'AlignmentIterator',
                              'AlignmentConstIterator',
                              'Scorer',
                              'Weightor',
                              'Renderer',
                              'LogOddor',
                              'Regularizor',
                              'Iterator2D',
                              'ResiduePAIR',
                              'Alignatum',
                              'EVDParameters',
                              'AlignandumData',
                              'SubstitutionMatrixData',
                              'NormalDistributionParameters',
                              'Distor',
                              'Treetor',
                              'Tree',
                              'PhyloMatrix',
                              ])

    ## include all classes
    mb.classes( lambda x: x.name in classes_to_export ).include()

    ## deal with copying/cloning members functions
    mb.member_functions( lambda x: x.name == "getClone" ).call_policies = \
                         return_value_policy( manage_new_object )
    mb.member_functions( lambda x: x.name == "getNew" ).call_policies = \
                         return_value_policy( manage_new_object )

    ## add __str__ function for functions having defined the '<<' operator
    ## This automatically maps std::ostream to a string
    mb.free_operators( lambda x: x.name == "operator<<" ).include()    
    
    ## used for dottor objects
    ## mb.member_functions( lambda x: x.name == "getPairs" ).call_policies = \
    ##                      return_value_policy( manage_new_object )
    # return reference to existing values - look at managing this better later
    mb.member_functions( lambda x: x.name == "getPointer" ).call_policies = \
                         return_value_policy( reference_existing_object )
    mb.member_functions( lambda x: x.name == "getReference" ).call_policies = \
                         return_value_policy( reference_existing_object )
    ## mb.member_functions( lambda x: x.name == "getRowIndices" ).call_policies = \
    ##                      return_value_policy( reference_existing_object )
    mb.member_functions( lambda x: x.name == "getRow" ).call_policies = \
                         return_value_policy( reference_existing_object )
    mb.member_functions( lambda x: x.name == "align" ).call_policies = \
                         return_value_policy( reference_existing_object )
    mb.member_functions( lambda x: x.name == "fragment" ).call_policies = \
                         return_value_policy( reference_existing_object )
    mb.member_functions( "decode" ).call_policies = \
                         return_value_policy( manage_new_object )
    mb.member_functions( "encode" ).call_policies = \
                         return_value_policy( manage_new_object )
    mb.member_functions( "calculateWeights" ).call_policies = \
                         return_value_policy( manage_new_object )
    mb.member_functions( "calculateMatrix" ).call_policies = \
                         return_value_policy( manage_new_object )                             
    mb.member_functions( "calculateTree" ).call_policies = \
                         return_value_policy( manage_new_object )                             
    mb.member_functions( lambda mem_fun: mem_fun.name.startswith( "getNodes" ) ).call_policies = \
                         return_value_policy( manage_new_object )


    ## exclude the following because of unhandled arguments/return types
    mb.member_functions( "fillProfile").exclude()
    mb.member_functions( "fillFrequencies").exclude()
    mb.member_functions( "getMatrix" ).exclude()

    ## get an error for this. For Dottor, now obsolete.
    ## mb.member_functions( lambda x: x.name in ("getRowIndices",) ).exclude()

    ## do not include the increment/decrement and dereference operators, because there is no equivalent in python
    ## exlude functions while testing. Need to map return types later.
    mb.member_operators( lambda x: x.name in ("operator++", "operator--", "operator*", "operator->", "operator()") ).exclude()

    ## do not export the internal iterator interfaces. This makes Alignment
    ## virtual and the wrapper will cause compilation to fail.
    mb.classes( lambda x: x.name in ("Iterator", "ConstIterator")).exclude()

    ## see http://article.gmane.org/gmane.comp.python.c++/10177
    ## this is used to remove the wrapper classes for purely virtual classes.
    for c in classes_to_export:
        cls = mb.class_( c )
        ## remove virtual member functions. This will result in a class to be
        ## not virtual and hence no wrapper.
        members = cls.decls( declarations.virtuality_type_matcher( declarations.VIRTUALITY_TYPES.PURE_VIRTUAL ),
                             decl_type=declarations.calldef.member_calldef_t,
                             allow_empty=True)
        members.set_virtuality( declarations.VIRTUALITY_TYPES.NOT_VIRTUAL )
        members = cls.decls( declarations.virtuality_type_matcher( declarations.VIRTUALITY_TYPES.ALL ),
                             decl_type=declarations.calldef.member_calldef_t,
                             allow_empty=True)
        members.set_virtuality( declarations.VIRTUALITY_TYPES.NOT_VIRTUAL )

        ## do not wrap constructors, because compilation will fail for
        ## abstract classes.
        cls.constructors().exclude()

    ## deal with methods that transfer ownership of the first argument
    ## supplied to the method.
    ## The list contains 
    ## 1: the class with the member function
    ## 2: the function prefix that needs to wrapped
    ## 3: the class of the pointee
    classes_with_ownership_transfer = [ ("MultipleAlignment", "add", "Alignatum") , ]
    ## TODO: add others, in particular addPair, but check for argument types 
                                         #("Alignment", "addPair", "ResiduePAIR" ) ]
    
    for ccontainer, fname, cpointee in classes_with_ownership_transfer:
        cls_pointee = mb.class_( cpointee )
        cls_pointee.held_type = 'std::auto_ptr< %s >' % cls_pointee.decl_string
                    
        cls_container= mb.class_( ccontainer )
        mem_funs = cls_container.member_functions( fname ) 
        for f in mem_funs:
            # set alias to original function name, otherwise ugly names will be created
            # a call for rename() had no effect.
            f.add_transformation( function_transformers.transfer_ownership( 0 ), 
                                  alias = fname )
                       
       
    ## export load/save functionality        
    exportSave( mb, options )        
    exportLoad( mb, options )
            
    ## Deal with templated matrix class
    template_translations = { 'Matrix<double>' : 'MatrixDouble',
                              'Matrix<unsigned>' : 'MatrixUInt',  
                              'Matrix<int>' : 'MatrixInt',  
                              }

    for old, new in template_translations.items():
        cls = mb.class_( old )
        cls.rename( new )
        cls.alias = new
        ## no warning for warning W1036: Py++ can not expose pointer to Python immutable member
        cls.vars(lambda x: x.name == "mMatrix" ).disable_warnings( messages.W1036 )

    declarations_to_export = set( template_translations.keys() )

    # mb.decls( lambda x: x.name in declarations_to_export ).include()


    
def buildModule( include_paths, dest, options) :
    """build module using py++."""
    
    if options.force:
        if os.path.exists( "cache" ):
            os.remove( "cache" )
            
    #Creating an instance of class that will help you to expose your declarations
    mb = module_builder.module_builder_t( [r"includes.h"]
                                          , gccxml_path=r""
                                          , cache="cache"
                                          , start_with_declarations=( "alignlib","py_details")
                                          , working_directory=r"."
                                          , include_paths=include_paths
                                          , define_symbols=[]
                                          , )
    
    ## Every declaration will be exposed at its own line
    mb.classes().always_expose_using_scope = True
    
    ## exclude py_details namespace, only used to instantiate template classes
    mb.namespace( 'py_details' ).exclude()

    addStreamBufClasses( mb )
        
    export_functions( mb )
    
    export_classes( mb )
    
    export_interface_classes( mb )
        
    ######################################################################
    
    # holder = mb.class_( 'vector<double>' )
    # holder.rename( 'VectorDouble' )
    
    
    ######################################################################
    #Well, don't you want to see what is going on?
    if options.verbose:
        mb.print_declarations()
    
    enumerations_to_export = set( ['AlignmentType', 'CombinationMode', 'SearchType', 'LinkageType' ] )
    
    mb.enumerations( lambda x: x.name in enumerations_to_export ).include()
    
    #Creating code creator. After this step you should not modify/customize declarations.
    mb.build_code_creator( module_name='alignlib' )
    mb.code_creator.add_include( "iostream" )
    mb.code_creator.add_include( "cstdio" )
    
    #Writing code to file.
    mb.write_module( dest )

def compileModule( module_name, options ):
    
    jamroot_template="""
# Copyright David Abrahams 2006. Distributed under the Boost
# Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

# Specify the path to the Boost project.  If you move this project,
# adjust this path to refer to the Boost root directory.
use-project boost 
  : %(boost_dir)s ;

# Set up the project-wide requirements that everything uses the
# boost_python library from the project whose global ID is
# /boost/python.
project
  : requirements 
    <library>/boost/python//boost_python 
    <library>%(alignlib_dir)s/libalignlib.so ;

# Declare a Python extension called hello.
python-extension alignlib 
  : # source 
    alignlib.cpp
  : # requirements 
  <library-file>%(alignlib_dir)s/libalignlib.so
  %(includes)s 
  ;
  """

    params = { "boost_dir" : options.boost_dir,
                "alignlib_dir" : options.alignlib_lib_dir,
                "includes" : "<include>" + "<include>".join( options.alignlib_include_dirs) }

    outfile = open("Jamroot", "w" )
    outfile.write(jamroot_template % params)
    outfile.close()
    
    s = subprocess.Popen( "bjam",
                    shell = True,
                    stdout = subprocess.PIPE,
                    stderr = subprocess.PIPE,
                    cwd = options.build_dir,
                    close_fds = True)

    (out, err) = s.communicate()

    print out
    print err
    
    if s.returncode != 0:
        raise "compilation error"

if __name__ == "__main__":
    
    parser = optparse.OptionParser( version = "%prog version: $Id$", usage = USAGE )

    parser.add_option( "-f", "--force", dest="force", action="store_true",
                      help="force complete rebuilt..")
    
    parser.add_option( "--boost-dir", dest="boost_dir", type="string",
                       help="location of boost." )
    
    parser.set_defaults( extension_name = "alignlib",
                         force = False, 
                         src_dir = "../src",
                         boost_dir = None,
                         alignlib_lib_dir = "../src/.libs",
                         alignlib_include_dirs = ["../src", ],
                         build_dir = ".",
                         verbose = False,
                         )
    
    (options, args) = parser.parse_args()

    if len(args) < 1:
        print USAGE
        raise "please supply a command"

    ## find boost
    findBoost( options )

    commands = map( lambda x: x.lower(), args)
    
    for command in commands:
        if command not in ("build", "test", "install", "generate-interface", "compile-interface" ):
            print USAGE
            raise "unknown command %s" % command
        
    for command in commands:

        if command in ("build", "generate-interface", "compile-interface" ):
    
            nerrors = checkRequisites( options )
            
            if nerrors:
                print "found %i errors - aborting build." % (nerrors)
            
            ## installation directory of alignlib
            src_dir=os.path.abspath( options.src_dir )
            
            module_name = "%s.cpp" % options.extension_name
            
            if options.force or not os.path.exists( module_name):
                print "building module %s" % module_name        
                buildModule( include_paths = [src_dir,], dest = module_name, options = options )
                
            if command == "generate-interface": break                
        
            print "compiling extension %s" % options.extension_name
        
            compileModule( module_name = module_name, options = options )

        elif command == "install":
            pass
        
        elif command == "test":
            pass
