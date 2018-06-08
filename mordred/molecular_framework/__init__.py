from rdkit import Chem
import networkx as nx


def get_molecular_framework(mol):
    rings = [set(s) for s in Chem.GetSymmSSSR(mol)]
    return _get_linkers(mol, rings), rings


def _get_linkers(mol, rings):
    G = nx.Graph()
    Rd = {i: ("R", Ri) for Ri, R in enumerate(rings) for i in R}
    R = list(set(Rd.values()))
    NR = len(R)

    for bond in mol.GetBonds():
        a = bond.GetBeginAtomIdx()
        b = bond.GetEndAtomIdx()
        a = Rd.get(a, ("A", a))
        b = Rd.get(b, ("A", b))

        G.add_edge(a, b)

    linkers = []
    for Ri, Rj in ((i, j) for i in range(NR) for j in range(i + 1, NR)):
        Ra, Rb = R[Ri], R[Rj]
        try:
            between = [n for n in nx.shortest_path(G, Ra, Rb) if n != Ra and n != Rb]
            if "R" in set(map(lambda v: v[0], between)):
                continue

            linkers.append({i for _, i in between})
        except nx.NetworkXNoPath:
            pass

    return linkers
