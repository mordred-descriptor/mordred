from __future__ import division

from collections import defaultdict

import numpy as np

from ._mesh import SphereMesh
from .._util import atoms_to_numpy
from .._atomic_property import GetAtomicNumber, vdw_radii


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

    def __init__(self, radiuses, xyzs, level=4):
        self.rads = radiuses
        self.rads2 = radiuses ** 2
        self.xyzs = xyzs
        self._gen_neighbor_list()
        self.sphere = SphereMesh(level).vertices.T

    def _gen_neighbor_list(self):
        r = self.rads[:, np.newaxis] + self.rads

        d = np.sqrt(np.sum((self.xyzs[:, np.newaxis] - self.xyzs) ** 2, axis=2))

        ns = defaultdict(list)
        for i, j in np.transpose(np.nonzero(d <= r)):
            if i == j:
                continue

            ns[i].append((j, d[i, j]))

        for _, l in ns.items():
            l.sort(key=lambda i: i[1])

        self.neighbors = ns

    def atomic_sa(self, i):
        r"""Calculate atomic surface area.

        :type i: int
        :param i: atom index

        :rtype: float
        """
        sa = 4.0 * np.pi * self.rads2[i]

        neighbors = self.neighbors.get(i)

        if neighbors is None:
            return sa

        XYZi = self.xyzs[i, np.newaxis].T

        sphere = self.sphere * self.rads[i] + XYZi
        N = sphere.shape[1]

        for j, _ in neighbors:
            XYZj = self.xyzs[j, np.newaxis].T

            d2 = (sphere - XYZj) ** 2
            mask = (d2[0] + d2[1] + d2[2]) > self.rads2[j]
            sphere = np.compress(mask, sphere, axis=1)

        return sa * sphere.shape[1] / N

    def surface_area(self):
        r"""Calculate all atomic surface area.

        :rtype: [float]
        """
        return [self.atomic_sa(i) for i in range(len(self.rads))]

    @classmethod
    def from_mol(cls, mol, conformer=-1, solvent_radius=1.4, level=4):
        r"""Construct SurfaceArea from rdkit Mol type.

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
        rs = atoms_to_numpy(lambda a: vdw_radii[a.GetAtomicNum()] + solvent_radius, mol)

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

        for atom in PDBParser().get_structure("", pdb).get_atoms():
            rs.append(vdw_radii[GetAtomicNumber(atom.element)] + solvent_radius)
            coords.append(atom.coord)

        return cls(np.array(rs), np.array(coords), level)
