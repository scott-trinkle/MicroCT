from setuptools import setup

setup(name='microct',
      version='0.1',
      description='Working with synchrotron microCT images from APS',
      url='https://github.com/scott-trinkle/MicroCT',
      author='Scott Trinkle',
      author_email='tscott.trinkle@gmail.com',
      license='MIT',
      packages=['microct'],
      install_requires=[
          'matplotlib',
          'numpy',
          'scipy',
          'libtiff',
          'netCDF4',
      ],
      zip_safe=False)
