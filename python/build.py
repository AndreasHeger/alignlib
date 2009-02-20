
USAGE = """python setup.py [options] command [command [...]]

The following commands are available:

* build: build and compile python extension to alignlib
* test: run some tests
* install: install the extension
"""

import re, sys, os, optparse, subprocess, glob

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

        for x in ( "/usr", "/usr/local/boost", "/opt/boost" ):
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
        print "could not find alignlib shared library %s/libalignlib.so" % options.alignlib_lib_dir
        
    return nerrors

def addStreamBufClasses( cls ):
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
    cls.add_declaration_code ( declaration_code )
    # mb.add_declaration_code ( declaration_code, tail = True )

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
    # mb.add_declaration_code( declaration_code, tail=True )

    for c in classes:
        cls = mb.class_(c)
        cls.include_files.append( "streambuf" )
        registration_code = 'def( "save", wrapper_for_save<alignlib::%s> )' % cls.name
        cls.member_function( "save" ).exclude()
        # the order is important
        addStreamBufClasses( cls )
        cls.add_declaration_code( declaration_code )
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
    for prefix in ("getMapResidue", "makeAlignatorDotsWrap", "getDefaultPalette" ):
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
                              'AlignmentFormatEmissions',
                              'MultAlignmentFormat',
                              'MultAlignmentFormatPlain' ,] )
    
    ## include all classes
    mb.classes( lambda x: x.name in classes_to_export ).include()

    ## add __str__ function for functions having defined the '<<' operator
    ## This automatically maps std::ostream to a string
    mb.free_operators( lambda x: x.name == "operator<<" ).include()    

    mb.member_functions( lambda x: x.name in ("applyOffset", "removeOffset" ) ).exclude()
    
    classes_with_load_save = set( ['AlignmentFormat',
                              'AlignmentFormatBlocks',
                              'AlignmentFormatExplicit',
                              'AlignmentFormatDiagonals',
                              'AlignmentFormatEmissions',
                              'MultAlignmentFormat',
                              'MultAlignmentFormatPlain', ] )

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
                              'MultAlignment',
                              'Alignator',
                              'MultipleAlignator',
                              'Iterator',
                              'Encoder',
                              'Fragmentor',
                              'Alignment',
                              'AlignmentIterator',
                              'Scorer',
                              'Weightor',
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
                              'Toolkit',
                              'Segment',
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
                         'HMultAlignment',
                         'HAlignator',
                         'HMultipleAlignator',
                         'HEncoder',
                         'HFragmentor',
                         'HAlignment',
                         'HScorer',
                         'HWeightor',
                         'HLogOddor',
                         'HRegularizor',
                         'HIterator2D',
                         'HAlignatum',
                         'HDistor',
                         'HTreetor',
                         'HTree',
                         'HPhyloMatrix',
                         'HFragmentVector',
                         'HToolkit',
                         'HSegmentVector', 
                         'HCountVector' ]

    #for c in mb.classes():
    #    print c.name
    #    print c.aliases

    for handle in handles_to_export:
        
        #for c in mb.classes( lambda x: handle[1:]in x.name ):
        #    print "class=", c.name 
        #for c in mb.decls( lambda x: handle[1:] in x.name):
        #    print "decl=", c.name
         
        # used to be: pointer = "boost::shared_ptr<alignlib::%s>" % handle[1:]
        pointer = "shared_ptr<alignlib::%s>" % handle[1:]
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



    handle_vectors_to_export = ( { 'name' : 'FragmentVector', 'handle' : 'HAlignment' }, 
                                 { 'name' : 'AlignandumVector', 'handle' : 'HAlignandum' } )

    for data in handle_vectors_to_export:
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


    atomic_vectors_to_export = ( { 'name' : 'StringVector', 'content' : 'std::string' },
                                 #{ 'name' : 'CountVector', 'content' : 'unsigned long' },
                                 { 'name' : 'CountVector', 'content' : 'alignlib::Count' },
#                                 { 'name' : 'NodeVector', 'content' : 'unsigned long' }, 
                                 ) 
     
    for data in atomic_vectors_to_export:
        code = """
        { //::std::vector<%(content)s, std::allocator<%(content)s> >
        typedef bp::class_< std::vector<%(content)s, std::allocator<%(content)s> > > %(name)s_exposer_t;
        %(name)s_exposer_t %(name)s_exposer = %(name)s_exposer_t( "%(name)s" );
        bp::scope %(name)s_scope( %(name)s_exposer );
        %(name)s_exposer.def( bp::vector_indexing_suite< ::std::vector<%(content)s, std::allocator<%(content)s> >, true >() );
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
        cls = mb.class_( old )
        cls.rename( new )
        cls.alias = new
        ## no warning for warning W1036: Py++ can not expose pointer to Python immutable member
        cls.vars(lambda x: x.name == "mMatrix" ).disable_warnings( messages.W1036 )
        ## do not wrap [], gives rise to "invalid application of 'sizeof' to incomplete type"        
        cls.member_operators( "operator[]" ).exclude()
        cls.mem_fun( "permuteRows" ).exclude()
        cls.mem_fun( "permuteCols" ).exclude()
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
                                   'LinkageType', 'AlphabetType', 'StorageType',
                                   'ToolkitType', 'AggregateType' ] )
    
    mb.enumerations( lambda x: x.name in enumerations_to_export ).include()
    
def buildModule( include_paths, dest, options) :
    """build module using py++."""
    
    if not GLOBAL_HAS_PYPLUSPLUS:
        raise ValueError( "can not build the interface, as py++ is not installed." )
    
    if options.force:
        if os.path.exists( "cache" ):
            os.remove( "cache" )
            
    #Creating an instance of class that will help you to expose your declarations
    mb = module_builder.module_builder_t( [r"includes.h"]
                                          , gccxml_path=r""
                                          , cflags = options.gccxml_options
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
    
    my_exception = mb.class_( 'AlignlibException' )
    my_exception.translate_exception_to_string( 'PyExc_RuntimeError', 'exc.what()')
    
    ######################################################################
    if options.verbose:
        print "# declarations after building interface."        
        mb.print_declarations()

    # creating code creator. After this step you should not modify/customize declarations.
    mb.build_code_creator( module_name='alignlib' )
    mb.code_creator.add_include( "iostream" )
    mb.code_creator.add_include( "cstdio" )
    mb.split_module( "modules" )

    # Writing code to file.
    mb.write_module( dest )

    # patch the output. Explict declarations were repeated
    # after splitting the module.
    
    lines = open( dest, "r").readlines()
    keep = True
    outfile = open( dest, "w" )
    pattern = "class std_obuf: public std::streambuf"
    
    x = 0
    for line in lines:
        x += 1
        if keep: outfile.write(line)

        if line.startswith( pattern ):
            if keep: keep = False
            else: break

    outfile.write( "".join(lines[x:]))        
    outfile.close()

if __name__ == "__main__":
    
    parser = optparse.OptionParser( version = "%prog version: $Id$", usage = USAGE )

    parser.add_option( "-f", "--force", dest="force", action="store_true",
                      help="force complete rebuilt..")
    
    parser.add_option( "--boost-dir", dest="boost_dir", type="string",
                       help="location of boost." )
    
    parser.add_option( "--verbose", dest="verbose", action="store_true",
                       help="output details." )
    
    parser.add_option( "--gccxml-options", dest="gccxml_options", type="string",
                        help="flags to be passed to gccxml." )
    
    parser.set_defaults( extension_name = "alignlib",
                         force = False, 
                         src_dir = "../alignlib",
                         boost_dir = None,
                         compiler = None,
                         gccxml_options = "",
                         alignlib_lib_dir = "../alignlib/.libs",
                         alignlib_include_dirs = [ "../alignlib", ],
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

        if command in ("build", "generate-interface" ):
    
            nerrors = checkRequisites( options )
            
            if nerrors:
                print "found %i errors - aborting build." % (nerrors)
            
            ## installation directory of alignlib
            src_dir=os.path.abspath( options.src_dir )
            
            module_name = "%s.cpp" % options.extension_name
            
            if options.force or not os.path.exists( module_name):
                print "building module %s" % module_name        
                buildModule( include_paths = [src_dir, 
                                              os.path.abspath( ".." ), 
                                              options.boost_dir ], 
                                              dest = module_name, 
                                              options = options )
                
            if command == "generate-interface": break                
        
            if command == "install":
            
                python_lib = distutils.sysconfig.get_python_lib()
                python_lib_data = python_lib + "/alignlib"
                results = []

                try:
                    if not os.path.exists( python_lib_data ):  
                        os.mkdir( python_lib_data ) 
                    shutil.copy( "exposed_decl.pypp.txt", python_lib_data )
                except IOError, msg:
                    print "installation failed: %s" % str(msg)
                    
        elif command == "test":
            pass

    