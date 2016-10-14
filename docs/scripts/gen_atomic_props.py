import load_path
from mordred import _atomic_property as prop

load_path.nop()


print('''
.. _atomic_properties:

atomic properties
=================
'''[1:])


for v in prop.getter_list:
    print(v.short)
    print(' ' * 4 + v.long)
