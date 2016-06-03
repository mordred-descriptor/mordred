from setuptools import setup, find_packages
import sys

install_requires = [
    'six>=1.10',
    'numpy>=1.10',
    'networkx>=1.10',
    'tqdm>=3.7.1',
    'click>=6.6',
]

if sys.version_info < (3, 4, 0):
    install_requires += ['enum34']

sandbox = {}
exec(open('mordred/_version.py').read(), sandbox, sandbox)

setup(
    name='mordred',
    version=sandbox['__version__'],
    packages=find_packages(),

    package_data={
        'mordred': ['data/*']
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
