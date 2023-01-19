from setuptools import setup
from Cython.Build import cythonize
from numpy.distutils.misc_util import Configuration

setup(
    name='planet integrator',
    ext_modules=cythonize("planet_new.pyx"),
    zip_safe=False,
)
config = Configuration('optimize')
config.add_extension('planet_new', sources=['planet_new.c'])