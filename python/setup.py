
USAGE = """python setup.py [options] command [command [...]]

The following commands are available:

* build: build and compile python extension to alignlib
* test: run some tests
* install: install the extension
"""

import re, sys, os, optparse, subprocess

from types import *

try:
    from pyplusplus import module_builder, messages, function_transformers
    from pyplusplus.module_builder import call_policies
    from pyplusplus.decl_wrappers import \
    return_value_policy, manage_new_object, copy_const_reference, reference_existing_object, return_self, return_arg

    from pygccxml import declarations
    GLOBAL_HAS_PYPLUSPLUS = True
except ImportError:
    GLOBAL_HAS_PYPLUSPLUS = False

import distutils.sysconfig
import os.path
import shutil

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

def exportSave( mb, classes, options, generic = True):
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

    ## add declaration code only once, as it is templated
    mb.add_declaration_code( declaration_code, tail=True )

    for c in classes:
        cls = mb.class_(c)
        cls.include_files.append( "streambuf" )
        registration_code = 'def( "save", wrapper_for_save<alignlib::%s> )' % cls.name
        cls.member_function( "save" ).exclude()
        cls.add_registration_code( registration_code )

def exportLoad( mb, classes, options, generic = True ):
    """export load functions.
    
    The load functions returns an empty smart pointer. This
    is then translated to the None object.
    
    If generic is true, the method name is appended to the
        function name (like loadAlignandum).
    """
    
    for c in classes:
        
        params = { 'class' : c }
        if generic:
            params['function'] = "load"
        else:
            params['function'] = "load%s" % c
        
        declaration_code = \
    """
      alignlib::H%(class)s wrapper_for_load_%(class)s( PyObject * fp )
      {
          if (!PyFile_Check(fp))
          {
              throw boost::python::error_already_set();
          }
         std::FILE * f = PyFile_AsFile(fp);   
         std_ibuf buf(f);
         std::istream is(&buf);
         if (is.peek() == EOF)
             return alignlib::H%(class)s();
         else
             return alignlib::H%(class)s (alignlib::%(function)s( is ));
      } 
""" % params

        registration_code = 'bp::def( "%(function)s", wrapper_for_load_%(class)s );' % params     
    
        fun = mb.free_function( "%(function)s" % params)
        fun.include_files.append( "streambuf" )
        fun.exclude()
        mb.add_declaration_code( declaration_code, tail = True )
        mb.add_registration_code( registration_code )

def exportFunctions( mb ):
    """export utility functions."""

    ## include all free functions
    mb.namespace("alignlib").free_functions().include()

    # change these functions to use handles.
    for prefix in ( "makeEVDParameters", "makeNormalDistributionParameters", "makeEntropyVector" ): 
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
   
def exportClasses( mb ):
    """export classes.
    
    These classes can be instantiated directly from python.
    """
    classes_to_export = set( ['Coordinate', 
                              'AlignmentFormat',
                              'AlignmentFormatBlocks',
                              'AlignmentFormatBlat',
                              'AlignmentFormatExplicit',
                              'AlignmentFormatDiagonals',
                              'AlignmentFormatEmissions', ] )
    
    ## include all classes
    mb.classes( lambda x: x.name in classes_to_export ).include()

    ## add __str__ function for functions having defined the '<<' operator
    ## This automatically maps std::ostream to a string
    mb.free_operators( lambda x: x.name == "operator<<" ).include()    
    
    classes_with_load_save = set( ['AlignmentFormat',
                              'AlignmentFormatBlocks',
                              'AlignmentFormatExplicit',
                              'AlignmentFormatDiagonals',
                              'AlignmentFormatEmissions', ] )

#    exportLoad( mb, 
#                classes_with_load_save,
#                options )
        
def exportInterfaceClasses( mb ):
    """export virtual classes.
    
    These classes can not instantiated directly from python.
    """
    classes_to_export = set( ['Alignandum',
                              'Sequence',
                              'Profile',
                              'MultipleAlignment',
                              'Alignator',
                              'Iterator',
                              'Encoder',
                              'Fragmentor',
                              'Alignment',
                              'AlignmentIterator',
                              'Scorer',
                              'Weightor',
                              'Renderer',
                              'LogOddor',
                              'Regularizor',
                              'Iterator2D',
                              'ResiduePair',
                              'Alignatum',
                              'EVDParameters',
                              'NormalDistributionParameters',
                              'Distor',
                              'Treetor',
                              'Tree',
                              'DistanceMatrix',
                              ])

    ## TODO export substitution matrix

    ## include all classes
    mb.classes( lambda x: x.name in classes_to_export ).include()

    ## do not include the increment/decrement and dereference operators, because there is no equivalent in python
    ## exlude functions while testing. Need to map return types later.
    ## default: exclude most operators
    mb.member_operators( lambda x: x.name in ("operator++", "operator--", 
                                              "operator*", "operator->", 
                                              "operator()", "operator[]") ).exclude()

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
        try:
            cls.constructors().exclude()
        except RuntimeError:
            print "no construtors for %s" % c
            

    ## export load/save functionality
    classes_with_load_save = ("Alignandum", "Encoder")         
    exportSave( mb, 
                classes_with_load_save,
                options,
                generic = False )        
    
    exportLoad( mb, 
                classes_with_load_save,
                options,
                generic = False )

def exportHandles( mb ):
    """include handle classes. 
    
    These are shared_ptr<> typedefs.     
    """
    handles_to_export = ['HAlignandum',
                         'HProfile',
                         'HSequence',
                         'HMultipleAlignment',
                        'HAlignator',
                        'HEncoder',
                        'HFragmentor',
                        'HAlignment',
                              'HScorer',
                              'HWeightor',
                              'HRenderer',
                              'HLogOddor',
                              'HRegularizor',
                              'HIterator2D',
                              'HAlignatum',
                              'HDistor',
                              'HTreetor',
                              'HTree',
                              'HPhyloMatrix',
                              'HFragmentVector',
                              ]

    for handle in handles_to_export:
        
        # for c in mb.classes( lambda x: handle[1:]in x.name ):
        #     print "class=", c.name 
        # for c in mb.decls( lambda x: handle[1:] in x.name):
        #    print "decl=", c.name
            
        pointer = "boost::shared_ptr<alignlib::%s>" % handle[1:]
        try:
            d = mb.decl( pointer )
        except RuntimeError:
            print "could not find handle class %s, searching with %s" % (handle, pointer)
            continue
        
        # print "exporting %s" % handle
        
        d.include()
        # d.rename( handle )
        # d.alias = handle

def exportContainers( mb ):
    """include containers of complex types. 
    
    Containers for atomic types worked out of the box, but
    I could not get the name mapping right for containers
    containing complex types like handle classes. I thus
    add explicit code.
    """

# The following did not work:
#    template_translations = { 'vector<boost::shared_ptr<alignlib::Alignment>, std::allocator<boost::shared_ptr<alignlib::Alignment> > >': 'FragmentVector',
#                             }
#      
#    ## first include the whole class, then exclude specific member functions
#    declarations_to_export = set( template_translations.keys() )
#    mb.decls( lambda x: x.name in declarations_to_export ).include()
#
#    for old, new in template_translations.items():
#        cls = mb.class_( old )
#        cls.rename( new )
#        cls.alias = new

    vectors_to_export = ( { 'name' : 'FragmentVector', 'handle' : 'HAlignment' }, )

    for data in vectors_to_export:
        
        code = """
        { //::std::vector<%(handle)s, std::allocator<%(handle)s> >
        typedef bp::class_< std::vector<alignlib::%(handle)s, std::allocator<alignlib::%(handle)s> > > %(name)s_exposer_t;
        %(name)s_exposer_t %(name)s_exposer = %(name)s_exposer_t( "%(name)s" );
        bp::scope %(name)s_scope( %(name)s_exposer );
        %(name)s_exposer.def( bp::vector_indexing_suite< ::std::vector<alignlib::%(handle)s, std::allocator<alignlib::%(handle)s> >, true >() );
        }

        bp::register_ptr_to_python< boost::shared_ptr<alignlib::%(name)s> >();
        """ % data

        mb.add_registration_code( code, tail=True )

def exportMatrices( mb ):
    """include matrix classes.
    """
    
    ## Deal with templated matrix class
    template_translations = { 'Matrix<double>' : 'MatrixDouble',
                              'Matrix<unsigned int>' : 'MatrixUInt',  
                              'Matrix<int>' : 'MatrixInt',
                              }

    exclude_functions = ( "getData", "setData", "copyData", "getRow" )
    
    ## first include the whole class, then exclude specific member functions
    declarations_to_export = set( template_translations.keys() )
    mb.decls( lambda x: x.name in declarations_to_export ).include()

    for old, new in template_translations.items():
        print old, new
        cls = mb.class_( old )
        cls.rename( new )
        cls.alias = new
        ## no warning for warning W1036: Py++ can not expose pointer to Python immutable member
        cls.vars(lambda x: x.name == "mMatrix" ).disable_warnings( messages.W1036 )
        ## do not wrap [], gives rise to "invalid application of 'sizeof' to incomplete type"        
        cls.member_operators( "operator[]" ).exclude()
        ## silence warnings about immutable types
        cls.mem_fun( "getValue" ).disable_warnings( messages.W1008 )
        cls.mem_fun( "setValue" ).disable_warnings( messages.W1008 )

        for f in exclude_functions:        
            try:
                cls.mem_fun( f ).exclude()
            except RuntimeError:
                continue

def exportEnums( mb ):
    """export enums."""
    enumerations_to_export = set( ['AlignmentType', 'CombinationMode', 'SearchType', 
                                   'LinkageType', 'AlphabetType', 'StorageType' ] )
    
    mb.enumerations( lambda x: x.name in enumerations_to_export ).include()
    
def buildModule( include_paths, dest, options) :
    """build module using py++."""
    
    if not GLOBAL_HAS_PYPLUSPLUS:
        raise "can not build the interface, as py++ is not installed."
    
    if options.force:
        if os.path.exists( "cache" ):
            os.remove( "cache" )
            
    #Creating an instance of class that will help you to expose your declarations
    mb = module_builder.module_builder_t( [r"includes.h"]
                                          , gccxml_path=r""
                                          , cache="cache"
                                          , start_with_declarations=( "alignlib","py_details" )
                                          , working_directory=r"."
                                          , include_paths=include_paths
                                          , define_symbols=[]
                                          , )
    
    ## exclude py_details namespace, only used to instantiate template classes
    mb.namespace( 'py_details' ).exclude()

    if options.verbose:
        print "# declarations before building interface."
        mb.print_declarations()

    addStreamBufClasses( mb )
    
    exportFunctions( mb )
    
    exportClasses( mb )
    
    exportInterfaceClasses( mb )

    exportHandles( mb )
    
    ## Every declaration will be exposed at its own line
    mb.classes().always_expose_using_scope = True

    exportContainers( mb )

    exportMatrices( mb )

    exportEnums( mb )
    
    ######################################################################
    #Well, don't you want to see what is going on?
    if options.verbose:
        print "# declarations after building interface."        
        mb.print_declarations()
    
        
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
    
    s = subprocess.Popen( "bjam release",
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
    
    parser.add_option( "--verbose", dest="verbose", action="store_true",
                       help="output details." )
    
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

        if command in ("build", "generate-interface", "compile-interface", "install" ):
    
            nerrors = checkRequisites( options )
            
            if nerrors:
                print "found %i errors - aborting build." % (nerrors)
            
            ## installation directory of alignlib
            src_dir=os.path.abspath( options.src_dir )
            
            module_name = "%s.cpp" % options.extension_name
            
            if options.force or not os.path.exists( module_name):
                print "building module %s" % module_name        
                buildModule( include_paths = [src_dir, options.boost_dir ], dest = module_name, options = options )
                
            if command == "generate-interface": break                
        
            print "compiling extension %s" % options.extension_name
        
            compileModule( module_name = module_name, options = options )

            if command == "install":
            
                python_lib = distutils.sysconfig.get_python_lib()
                    
                results = []
                for root, dirs, files in os.walk('.'):
        
                    if "alignlib.so" in files:
                        results.append( os.path.join(root, "alignlib.so") )
                        
                if len(results) == 0:
                    print "could not find alignlib.so"
                    print "please run setup.py build first."
                    sys.exit(0)
                
                if len(results) > 1:
                    print "found more than one alignlib.so in %s." % str(results)
                    print "confused and thus not installing."
            
                print "installing %s in %s" % (results[0], python_lib )
                try:
                    shutil.copy( results[0], python_lib )
                except IOError, msg:
                    print "installation failed: %s" % str(msg)
                    
        elif command == "test":
            pass
