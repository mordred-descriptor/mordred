"""Calculate molecular framework."""

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
        i = bond.GetBeginAtomIdx()
        j = bond.GetEndAtomIdx()
        a = Rd.get(i, ("A", i))
        b = Rd.get(j, ("A", j))

        G.add_edge(a, b, bond=(i, j))

    linkers = []
    for Ri, Rj in ((i, j) for i in range(NR) for j in range(i + 1, NR)):
        Ra, Rb = R[Ri], R[Rj]
        try:
            path = list(nx.shortest_path(G, Ra, Rb))
            if "R" in {n[0] for n in path if n != Ra and n != Rb}:
                continue

            bonds = []
            for a, b in zip(path, path[1:]):
                bonds.append(G[a][b]["bond"])

            linkers.append(bonds)
        except nx.NetworkXNoPath:
            pass

    return linkers
