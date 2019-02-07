# setup.py is like Python makefile
# from Cython Basics, a clearer idea can be formed:
# http://cython.readthedocs.io/en/latest/src/tutorial/cython_tutorial.html

from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import numpy

ext_modules = [
	Extension("WAM_functions",["WAM_functions.pyx"],
	extra_compile_args=['-fopenmp'],
	extra_link_args=['-fopenmp'],
	)
]

setup(
    ext_modules=cythonize(ext_modules),
    include_dirs=[numpy.get_include()]
)
# include_dirs - to specify include directories to search;
# numpy.get_include() - return the directory that contains the NumPy \*.h header files. Extension modules that need to compile against NumPy should use this function to locate the appropriate include directory.
