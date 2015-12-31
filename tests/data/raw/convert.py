import sys
import csv
import yaml
try:
    from yaml import CLoader as Loader
except ImportError:
    from yaml import Loader


def parse_int_or_float(s):
    if not isinstance(s, str):
        return s

    try:
        return int(s)
    except ValueError:
        return float(s)


def main(file_base):
    smi_file = file_base + '.smi'
    rename_file = file_base + '.rename'
    override_file = file_base + '.yaml'
    csv_file = file_base + '.csv'

    # read smiles file
    name_smi = dict()
    for line in open(smi_file):
        words = line.strip().split()
        name_smi[' '.join(words[1:])] = words[0]

    # read rename file
    rename = dict()
    for line in open(rename_file):
        src, dst = line.strip().split()
        rename[src] = dst

    # read override file
    overrides = yaml.load(open(override_file), Loader=Loader)

    # read csv file
    reader = csv.reader(open(csv_file))
    header = next(reader)[1:]

    for row in reader:
        name = row[0]
        smi = name_smi[name]

        override = overrides.get(name, {})

        print('- name: {}'.format(name))
        print('  smiles: {}'.format(smi))
        print('  descriptors:')

        for name, value in zip(header, row[1:]):
            name = rename.get(name, name)
            value = override.get(name, value)

            print('    {}: {!r}'.format(name, parse_int_or_float(value)))

        print()


main(sys.argv[1])
