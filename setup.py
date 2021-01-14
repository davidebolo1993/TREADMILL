from setuptools import setup,find_packages


setup(name='TREADMILL',      
  version='1.0.0',
  description='Tandem REpeats AnD MethylatIon caLLing',
  url='https://github.com/davidebolo1993/TREADMILL',
  requires=['python (>= 3.8)'],
  author='Davide Bolognini',
  author_email='davidebolognini7@gmail.com',
  license='LICENSE.txt',
  install_requires=['pyfaidx>=0.5.5.2', 'pysam>=0.15.0', 'editdistance>=0.5.2', 'numpy>=1.15.3', 'pybedtools>=0.8.0', 'mappy>=2.17', 'pandas>=1.1.2', 'weighted-levenshtein>=0.2.1', 'scikit-learn>=0.24.0'],
  zip_safe=False,
  packages=find_packages(),
  include_package_data=True,
  entry_points={'console_scripts': ['TREADMILL=TREADMILL.TREADMILL:main']})
