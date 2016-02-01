import os
from rdkit import Chem
from mordred._surface_area.sasa import SurfaceArea
from nose.tools import ok_


# calculated by pymol
txt_data = '''
Hexane               296.910
Benzene              243.552
Caffeine             369.973
Cyanidin             483.873
Lycopene            1172.253
Epicatechin          489.498
Limonene             361.278
Allicin              356.872
Glutathione          530.679
Digoxin             1074.428
Capsaicin            641.527
EllagicAcid          440.267
Astaxanthin         1080.941
DMSO                 227.926
DiethylThioketone    290.503
VinylsulfonicAcid    246.033
Thiophene            227.046
Triethoxyphosphine   396.482
MethylphosphonicAcid 235.685
MethylCyclopropane   229.071
Acetonitrile         182.197
Histidine            335.672
'''[1:-1]


sdf_file = os.path.join(
    os.path.dirname(__file__),
    'references',
    'structures.sdf',
)


def test_SASA():
    data = {}
    for line in txt_data.split('\n'):
        n, v = line.strip().split()
        data[n] = float(v)

    for mol in Chem.SDMolSupplier(sdf_file, removeHs=False):
        name = mol.GetProp('_Name')
        actual = sum(SurfaceArea.from_mol(mol).surface_area())
        desired = data[name]

        e = actual / desired
        p = 0.05

        yield ok_, 1 - p < e < 1 + p, 'large SASA error in {}: {}'.format(name, e)

'''
if __name__ == '__main__':
    from rdkit import Chem

    import sys

    f = sys.argv[1]

    for mol in Chem.SDMolSupplier(f, removeHs=False):
        sasa = SurfaceArea.from_mol(mol, level=5).surface_area()
        print(mol.GetProp('_Name'), np.sum(sasa))
        '''
