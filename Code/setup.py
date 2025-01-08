from distutils.core import setup, Extension
from Cython.Build import cythonize
import numpy


setup(name="icl_work", ext_modules=cythonize('All_data.pyx'),include_dirs=[numpy.get_include()])

