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
    install_requires.append('enum34')


def get_version():
    version_file = os.path.join('mordred', '_version.py')

    sandbox = {}
    exec(open(version_file).read(), sandbox, sandbox)
    return sandbox['__version__']


def get_test_data():
    for p, _, fs in os.walk(os.path.join('mordred', 'tests', 'references')):
        p = p.split(os.sep)[2:]

        for f in fs:
            yield os.path.join(*(p + [f]))

setup(
    name='mordred',
    version=get_version(),
    packages=find_packages(),

    package_data={
        'mordred': ['data/*.txt'],
        'mordred.tests': list(get_test_data()),
    },

    install_requires=install_requires,

    tests_require=[
        'nose>=1.3',
        'PyYaml>=3.11',
    ],

    test_suite='nose.collector',

    author='Hirotomo Moriwaki',
    author_email='philopon.dependence@gmail.com',
    license='BSD3',
)
