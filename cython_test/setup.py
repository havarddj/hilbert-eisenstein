#!/home/havard/sage/sage/local/bin/python
from distutils.core import setup
from Cython.Build import cythonize

setup(
    name = 'Cython test library',
    # ext_modules = cythonize( 'library.pyx' ),
    ext_modules =("RM-set.spyx")
)
