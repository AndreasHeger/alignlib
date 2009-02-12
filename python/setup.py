from distutils.core import setup
from distutils.extension import Extension
import os.path
import sys, glob

include_dirs = ["/usr/include",".","../alignlib", ".."]
libraries=["boost_python"]
library_dirs=['/usr/lib']

files = glob.glob( "modules/*.cpp" )

setup(name="alignlib",    
      ext_modules=[
                   Extension("alignlib",
                             files,
                             library_dirs=library_dirs,
                             libraries=libraries,
                             include_dirs=include_dirs,
                             depends=[]),
                             ]
     )
