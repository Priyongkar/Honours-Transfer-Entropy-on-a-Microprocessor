from setuptools import setup, Extension
from Cython.Build import cythonize
import numpy as np
ext = Extension('sim2600.mycircuitsimulator', ['sim2600/mycircuitsimulator.pyx', 'sim2600/cirsim.cc'],
                language='c++', include_dirs=[np.get_include(), 'sim2600'], extra_compile_args=['-O2', '-std=c++11'])
setup(name='sim2600', packages=['sim2600'], ext_modules=cythonize([ext], language_level=3))
