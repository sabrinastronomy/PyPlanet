from distutils.core import setup
from Cython.Build import cythonize

setup(ext_modules = cythonize('Pre_calc_liquidpv_mantle_density.pyx'))
