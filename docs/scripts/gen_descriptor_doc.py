import sys
from itertools import chain

import load_path
from mordred import descriptors
from mordred._base.descriptor import is_descriptor_class

load_path.nop()


class Generator(object):
    def __init__(self, mdl):
        self.mdl = mdl

    def header(self):
        h = self.mdl.__name__ + ' module'
        yield h
        yield '=' * len(h)
        yield ''

    def descriptor(self, desc):
        excludes = set([
            'as_argument', 'calculate', 'coord', 'dependencies', 'fail',
            'is_descriptor_class', 'mol', 'parameters', 'preset',
            'rethrow_na', 'rethrow_zerodiv'
        ])
        docs = [m for m in dir(desc) if m[:1] != '_' if m not in excludes]

        yield '.. autoclass:: {}.{}'.format(self.mdl.__name__, desc.__name__)
        yield '    :members: {}'.format(', '.join(docs))
        yield '    :undoc-members:'
        yield '    :show-inheritance:'
        yield ''

    def member(self, member):
        if is_descriptor_class(member):
            for line in self.descriptor(member):
                yield line
        else:
            raise ValueError('unknown member {} in {}'.format(member, self.mdl.__name__))

    def members(self):
        for m in self.mdl.__all__:
            for line in self.member(getattr(self.mdl, m)):
                yield line

    def module(self):
        yield '.. automodule:: {}'.format(self.mdl.__name__)

    def print(self, file=sys.stdout):
        for line in chain(self.header(), self.module(), self.members()):
            print(line, file=file)


if __name__ == '__main__':
    print('Submodules')
    print('----------')
    print('')
    print('.. toctree::')

    for mdl in descriptors.all:
        print('    ' + mdl.__name__)
        with open('api/{}.rst'.format(mdl.__name__), 'w') as f:
            Generator(mdl).print(file=f)
