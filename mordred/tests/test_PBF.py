# Based on the code from https://github.com/rdkit/rdkit/blob/master/Contrib/PBF/pbf.py
import numpy as np
from numpy import linalg
from rdkit import Chem


def test_PBF():
    def GetBestFitPlane(pts, weights=None):
        if weights is None:
            wSum = len(pts)
            origin = np.sum(pts, 0)
        origin /= wSum

        sums = np.zeros((3, 3), np.double)
        for pt in pts:
            dp = pt - origin
            for i in range(3):
                sums[i, i] += dp[i] * dp[i]
                for j in range(i + 1, 3):
                    sums[i, j] += dp[i] * dp[j]
                    sums[j, i] += dp[i] * dp[j]
        sums /= wSum
        vals, vects = linalg.eigh(sums)
        order = np.argsort(vals)
        normal = vects[:, order[0]]
        plane = np.zeros((4,), np.double)
        plane[:3] = normal
        plane[3] = -1 * normal.dot(origin)
        return plane

    def PBFRD(mol, confId=-1):
        conf = mol.GetConformer(confId)
        if not conf.Is3D():
            return 0

        pts = np.array([list(conf.GetAtomPosition(x)) for x in range(mol.GetNumAtoms())])
        plane = GetBestFitPlane(pts)
        denom = np.dot(plane[:3], plane[:3])
        denom = denom ** 0.5
        # add up the distance from the plane for each point:
        res = 0.0
        for pt in pts:
            res += np.abs(pt.dot(plane[:3]) + plane[3])
        res /= denom
        res /= len(pts)
        return res

    suppl = Chem.SDMolSupplier('./testData/egfr.sdf', removeHs=False)
    expected = open('./testData/egfr.out', 'r')
    for m in suppl:
        res = PBFRD(m)
        inl = next(expected).strip().split()
        expect = float(inl[1])
        assert abs(res - expect) < 1e-4
