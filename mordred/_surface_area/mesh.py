import math

import numpy as np

# http://prideout.net/blog/?p=44
# http://blog.andreaskahler.com/2009/06/creating-icosphere-mesh-in-code.html


class Vector3D(object):
    def __init__(self, x, y, z):
        self.x, self.y, self.z = x, y, z

    def __add__(self, other):
        return self.__class__(
            self.x + other.x,
            self.y + other.y,
            self.z + other.z,
        )

    def __abs__(self):
        return math.sqrt(self.x ** 2 + self.y ** 2 + self.z ** 2)

    def normalize(self):
        d = abs(self)

        if d:
            self.x /= d
            self.y /= d
            self.z /= d

        return self

    def __iter__(self):
        return iter([self.x, self.y, self.z])


class SphereMesh(object):
    def __init__(self, level=4):
        t = (1.0 + math.sqrt(5.0)) / 2.0
        self.vertices = [
            Vector3D(-1, t, 0).normalize(),
            Vector3D(1, t, 0).normalize(),
            Vector3D(-1, -t, 0).normalize(),
            Vector3D(1, -t, 0).normalize(),

            Vector3D(0, -1, t).normalize(),
            Vector3D(0, 1, t).normalize(),
            Vector3D(0, -1, -t).normalize(),
            Vector3D(0, 1, -t).normalize(),

            Vector3D(t, 0, -1).normalize(),
            Vector3D(t, 0, 1).normalize(),
            Vector3D(-t, 0, -1).normalize(),
            Vector3D(-t, 0, 1).normalize(),
        ]

        self.faces = [
            (0, 11, 5),
            (0, 5, 1),
            (0, 1, 7),
            (0, 7, 10),
            (0, 10, 11),

            (1, 5, 9),
            (5, 11, 4),
            (11, 10, 2),
            (10, 7, 6),
            (7, 1, 8),

            (3, 9, 4),
            (3, 4, 2),
            (3, 2, 6),
            (3, 6, 8),
            (3, 8, 9),

            (4, 9, 5),
            (2, 4, 11),
            (6, 2, 10),
            (8, 6, 7),
            (9, 8, 1),
        ]

        self.level = 1
        self.subdivide(level)

    def _subdivide(self):
        self.level += 1
        vs = self.vertices
        fs = self.faces

        for faceIndex, face in list(enumerate(fs)):
            a, b, c = (vs[i] for i in face)

            vs.append((a + b).normalize())
            vs.append((b + c).normalize())
            vs.append((a + c).normalize())

            i = len(vs) - 3
            j, k = i + 1, i + 2

            fs.append((i, j, k))
            fs.append((face[0], i, k))
            fs.append((i, face[1], j))
            fs[faceIndex] = (k, j, face[2])

    def subdivide(self, level):
        for _ in range(level - self.level):
            self._subdivide()

    @property
    def vertices_numpy(self):
        return np.array([list(v) for v in self.vertices], dtype='float')

    def __len__(self):
        return len(self.vertices)
