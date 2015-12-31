from setuptools import setup, find_packages
import sys

install_requires = [
    'six>=1.10',
    'numpy>=1.10',
    'scipy>=0.16',
    'networkx>=1.10',
]

if sys.version_info < (3, 4, 0):
    install_requires += ['enum34']

setup(
    name='mold',
    version='0.0.0',
    packages=find_packages(),

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
