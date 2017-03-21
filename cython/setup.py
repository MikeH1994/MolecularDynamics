from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

ext_modules = [Extension('pymoldC', ['pymold.pyx'],)]

setup(
    name="Set 1 of Functions",
    cmdclass={'build_ext': build_ext},
    ext_modules=ext_modules
)
