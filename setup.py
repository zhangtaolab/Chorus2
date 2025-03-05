# -*- coding: utf-8 -*-

import os
from setuptools import setup
from setuptools import find_packages

def readme(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(name='Chorus2',
      version='2.0',
      packages=find_packages(),
      author='Tao Zhang',
      author_email='zhangtao@yzu.edu.cn',
      scripts=['bin/Chorus2', 'bin/ChorusNGSfilter', 'bin/ChorusNGSselect', 'bin/ChorusHomo',
               'bin/ChorusNoRef', 'bin/ChorusGUI', 'bin/ChorusPBGUI', 'bin/ChorusDraftPrebuild'],
      py_modules=['Chorus2', 'ChorusNGSfilter', 'ChorusNGSselect', 'ChorusHomo',
                  'ChorusNoRef', 'ChorusDraftPrebuild', 'ChorusGUI', 'ChorusPBGUI'],
      include_package_data=True,
      url='https://github.com/zhangtaolab/Chorus2',
      license='LICENSE',
      description='Chorus2 software for oligo-FISH probe design',
      long_description=readme('README.md'),
      classifiers=['Intended Audience :: Science/Research'],
      install_requires=["cython", "numpy", "pyfasta", "primer3-py", "pandas", "pyBigWig",
                        "pybedtools", "matplotlib", "PyQt5<5.11", "SIP>=4"],
      zip_safe=False,
     )
