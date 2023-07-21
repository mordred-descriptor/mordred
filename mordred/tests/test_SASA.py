import os

from rdkit import Chem

from mordred.CPSA import TotalSurfaceArea
from mordred.surface_area import SurfaceArea

# calculated by pymol
txt_data = """
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
"""[
    1:-1
]


sdf_file = os.path.join(os.path.dirname(__file__), "references", "structures.sdf")


def test_SASA():
    data = {}
    for line in txt_data.split("\n"):
        n, v = line.strip().split()
        data[n] = float(v)

    tsa = TotalSurfaceArea()

    for mol in Chem.SDMolSupplier(sdf_file, removeHs=False):
        name = mol.GetProp("_Name")
        actual = sum(SurfaceArea.from_mol(mol).surface_area())

        if name not in data:
            continue

        desired = data[name]

        e = abs((actual - desired) / desired)
        p = 0.05

        assert e < p, "large SASA error in {}: {}".format(name, e)

        e = abs((tsa(mol) - desired) / desired)

        assert e < p, "large SASA error in {}: {}".format(name, e)
