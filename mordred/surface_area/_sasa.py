from collections import defaultdict

import numpy as np

from ._mesh import SphereMesh
from .._atomic_property import table, Rvdw
from .._util import atoms_to_numpy


class SurfaceArea(object):
    r"""calculate solvent accessible surface area.

    :type radiuses: np.ndarray(dtype=float, shape=(N,))
    :param radiuses: atomic radius + solvent radius vector

    :type xyzs: np.ndarray(dtype=float, shape=(N, 3))
    :param xyzs: atomic position matrix

    :type level: int
    :param level: mesh level. subdivide icosahedron n-1 times.

        .. math::

            N_{\rm points} = 5 \times 4^{level} - 8
    """

    def __init__(self, radiuses, xyzs, level=5):
        self.rads = radiuses
        self.xyzs = xyzs
        self._gen_neighbor_list()
        self.sphere = SphereMesh(level=level).vertices_numpy

    def _gen_neighbor_list(self):
        r = self.rads[:, np.newaxis] + self.rads

        d = np.sqrt(np.sum(
            (self.xyzs[:, np.newaxis] - self.xyzs) ** 2,
            axis=2))

        ns = defaultdict(list)
        for i, j in np.transpose(np.nonzero(d <= r)):
            if i == j:
                continue

            ns[i].append((j, d[i, j]))

        for _, l in ns.items():
            l.sort(key=lambda i: i[1])

        self.neighbors = ns

    def atomic_sa(self, i):
        r"""calculate atomic surface area.

        :type i: int
        :param i: atom index

        :rtype: float
        """

        Ri = self.rads[i]
        sa = 4.0 * np.pi * Ri ** 2

        neighbors = self.neighbors.get(i)

        if neighbors is None:
            return sa

        XYZi = self.xyzs[i]

        sphere = self.sphere * Ri + XYZi
        N = len(sphere)

        for j, _ in neighbors:
            Rj = self.rads[j]
            XYZj = self.xyzs[j]
            sphere = sphere[np.sqrt(np.sum((sphere - XYZj) ** 2, axis=1)) > Rj]

        return sa * float(len(sphere)) / N

    def surface_area(self):
        r"""calculate all atomic surface area.

        :rtype: [float]
        """

        return [self.atomic_sa(i) for i in range(len(self.rads))]

    @classmethod
    def from_mol(cls, mol, conformer=-1, solvent_radius=1.4, level=5):
        r"""construct SurfaceArea from rdkit Mol type.

        :type mol: rdkit.Chem.Mol
        :param mol: input molecule

        :type conformer: int
        :param conformer: conformer id

        :type solvent_radius: float
        :param solvent_radius: solvent radius

        :type level: int
        :param level: mesh level

        :rtype: SurfaceArea
        """

        rs = atoms_to_numpy(lambda a: Rvdw[a.GetAtomicNum()] + solvent_radius, mol)

        conf = mol.GetConformer(conformer)

        ps = np.array([list(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())])

        return cls(rs, ps, level)

    @classmethod
    def from_pdb(cls, pdb, solvent_radius=1.4, level=3):
        try:
            from Bio.PDB import PDBParser
        except ImportError:
            raise ImportError("There isn't biopython package.")

        rs = []
        coords = []

        for atom in PDBParser().get_structure('', pdb).get_atoms():
            rs.append(Rvdw[table.GetAtomicNumber(atom.element)] + solvent_radius)
            coords.append(atom.coord)

        return cls(np.array(rs), np.array(coords), level)
