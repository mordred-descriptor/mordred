import conf
from mordred import _atomic_property as prop


print('''
.. _atomic_properties:

atomic properties
=================
'''[1:])


for v in prop.getter_list:
    print(v.short)
    print(' ' * 4 + v.long)
