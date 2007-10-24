import re, sys, os
from pyplusplus import module_builder, messages
from pyplusplus.module_builder import call_policies
from pyplusplus.decl_wrappers import \
     return_value_policy, manage_new_object, copy_const_reference, reference_existing_object, \
     return_self, return_arg

from pygccxml import declarations

## installation directory of alignlib
src_dir="/home/andreas/devel"

#Creating an instance of class that will help you to expose your declarations
mb = module_builder.module_builder_t( [r"includes.h"]
                                      , gccxml_path=r""
                                      , cache="cache"
                                      , start_with_declarations=("alignlib","py_details")
                                      , working_directory=r"."
                                      , include_paths=['%s/include' % src_dir, '%s/include/alignlib' % src_dir]
                                      , define_symbols=[]
                                      , )

## Every declaration will be exposed at its own line
mb.classes().always_expose_using_scope = True

## exclude py_details namespace, only used to instantiate template classes
mb.namespace( 'py_details' ).exclude()

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
                   "perform",
                   "copy",
                   "combine",
                   "complement",
                   "rescore",
                   "flatten",
                   "filter",
                   "readAlignataPairs",
                   "calculateAffineScore",
                   "extractMultipleAlignment",
                   "rescaleProfileCounts",
                   "normalizeProfileCounts",
                   ):
        mb.free_functions( lambda mem_fun: mem_fun.name.startswith( prefix )).call_policies = return_self()

    # set call policies for functions that return get a default object
    # in this case the caller is not a new object.
    for prefix in ("getDefault", ):
        mb.free_functions( lambda mem_fun: mem_fun.name.startswith( prefix )).call_policies = \
                           return_value_policy( reference_existing_object )

    # set call policies for functions that return set a default object
    # in this case the caller takes ownership of the old object
    for prefix in ("setDefault", ):
        mb.free_functions( lambda mem_fun: mem_fun.name.startswith( prefix )).call_policies = \
                           return_value_policy( manage_new_object )

    # other functions that return new objects, including the factory functions
    for prefix in ("extract", "read", "make", "exportProfileFrequencies", "splitAlignata" ):
        mb.free_functions( lambda mem_fun: mem_fun.name.startswith( prefix )).call_policies = \
                           return_value_policy( manage_new_object )
        
    #######################################################################################
    #######################################################################################
    #######################################################################################
    ## patches to exclude problematic functions
    for prefix in ("getMapResidue", "makeSubstitutionMatrixAA", "makeRendererColumn", "makeAlignatorDotsWrap", "getDefaultPalette" ):
        mb.free_functions( lambda mem_fun: mem_fun.name.startswith( prefix )).exclude()

    ## can't export makeProfile(const std::string &, int)
    ## couldn't figure out what to mach const std::string & with
    ## tried: const std::string &, std::string const &, and more
    mb.free_functions( name='makeProfile', arg_types=[None, "int"] ).exclude()
    
        
def export_classes( mb ):
    """export classes."""
    classes_to_export = set( ['Alignandum',
                              'MultipleAlignment',
                              'Alignator',
                              'Iterator',
    #                          'Dottor',
    #                          'Translator',
                              'SubstitutionMatrix',
                              'Fragmentor',
                              'Alignata',
                              'AlignataIterator',
                              'AlignataConstIterator',
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

    ## exclude the following because of unhandled arguments/return types
    mb.member_functions( "fillProfile").exclude()
    mb.member_functions( "fillFrequencies").exclude()
    mb.member_functions( "getMatrix" ).exclude()

    ## get an error for this. For Dottor, now obsolete.
    ## mb.member_functions( lambda x: x.name in ("getRowIndices",) ).exclude()

    ## do not include the increment/decrement and dereference operators, because there is no equivalent in python
    ## exlude functions while testing. Need to map return types later.
    mb.member_operators( lambda x: x.name in ("operator++", "operator--", "operator*", "operator->", "operator()") ).exclude()

    ## do not export the internal iterator interfaces. This makes Alignata
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

    ## Deal with templated matrix class
    template_translations = { 'Matrix<double>' : 'MatrixDouble',
                              'Matrix<unsigned>' : 'MatrixUInt',  
                              'Matrix<int>' : 'MatrixInt',  
                              }

    for old, new in template_translations.items():
        cls = mb.class_( old )
        cls.rename( new )
        ## no warning for warning W1036: Py++ can not expose pointer to Python immutable member
        cls.vars(lambda x: x.name == "mMatrix" ).disable_warnings( messages.W1036 )

    declarations_to_export = set( template_translations.keys() )

    # mb.decls( lambda x: x.name in declarations_to_export ).include()

export_functions( mb )

export_classes( mb )

######################################################################

# holder = mb.class_( 'vector<double>' )
# holder.rename( 'VectorDouble' )

######################################################################
#Well, don't you want to see what is going on?
mb.print_declarations()

enumerations_to_export = set( ['AlignmentType', 'CombinationMode'] )

mb.enumerations( lambda x: x.name in enumerations_to_export ).include()

#Creating code creator. After this step you should not modify/customize declarations.
mb.build_code_creator( module_name='alignlib' )

#Writing code to file.
mb.write_module( './alignlib.cpp' )
