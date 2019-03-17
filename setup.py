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
      scripts=['bin/Chorus2', 'bin/ChorusNGSfilter', 
               'bin/ChorusNGSselect', 'bin/ChorusDraftPrebuild'],
      py_modules=['Chorus2', 'ChorusNGSfilter', 'ChorusNGSselect', 'ChorusDraftPrebuild'],
      include_package_data=True, 
      url='https://github.com/zhangtaolab/Chorus2',
      description='Chorus2 Software for Oligo FISH probe design',
      long_description=readme('README.md'),
      install_requires=["cython", "numpy", "pyfasta", "primer3-py", "pandas", "pyBigWig"],
      zip_safe=False,
      platforms=['Linux', 'OS X']
     )
