from setuptools import setup, find_packages
import sys
import os

install_requires = [
    'six>=1.10',
    'numpy>=1.10',
    'networkx>=1.10',
    'tqdm>=3.7.1',
]

if sys.version_info < (3, 4, 0):
    install_requires += ['enum34']

setup(
    name='mordred',
    version=open(os.path.join(os.path.dirname(__file__), 'mordred', 'version.txt')).read().strip(),
    packages=find_packages(),

    package_data={
        'mordred': ['version.txt']
    },

    install_requires=install_requires,

    extras_require=dict(
        test=['nose>=1.3',
              'PyYaml>=3.11',
              ]
    ),

    test_suite='nose.collector',

    author='Hirotomo Moriwaki',
    author_email='philopon.dependence@gmail.com',
    license='BSD3',
)
