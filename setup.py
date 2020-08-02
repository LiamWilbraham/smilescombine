from distutils.core import setup
from setuptools import find_packages

setup(
  name = 'smilescombine',
  packages = ['smilescombine'],
  version = '0.2',
  license='MIT',
  description = 'Combinatorially combine molecular skeletons and substituents.',
  author = 'Liam Wilbraham',
  author_email = 'liam.wilbrahaml@glasgow.ac.uk',
  url = 'https://github.com/LiamWilbraham/smilescombine',
  download_url = 'https://github.com/LiamWilbraham/smilescombine/archive/v_02.tar.gz',
  keywords = ['cheminformatics', 'chemistry'],
  classifiers=[
    'License :: OSI Approved :: MIT License',
    'Programming Language :: Python :: 3',
  ],
)