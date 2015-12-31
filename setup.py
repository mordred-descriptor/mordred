from setuptools import setup, find_packages
import sys

install_requires = [
    'six',
    'numpy',
    'scipy',
    'networkx',
]

if sys.version_info < (3, 4, 0):
    install_requires += ['enum34']

setup(
    name='mold',
    version='0.0.0',
    packages=find_packages(),

    install_requires=install_requires,

    extras_require=dict(
        test=['nose', 'PyYaml']
    ),

    test_suite='nose.collector',

    author='Hirotomo Moriwaki',
    author_email='philopon.dependence@gmail.com',
    license='BSD3',
)
