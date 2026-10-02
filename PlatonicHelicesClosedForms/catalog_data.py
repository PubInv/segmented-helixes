"""Inputs transcribed from Read (2022), Table 1; unit-edge solid geometry.

Standard library only. Face numbering follows repository revision
8be3839de2204b966c8aedb7c28f6a6f5c18076e. The companion verifier distinguishes
the printed catalog twists from analytically corrected vertex-matching twists.
"""
from itertools import product
from math import sqrt

SOURCE_REVISION = "8be3839de2204b966c8aedb7c28f6a6f5c18076e"
# name, solid, analogs, C face, twist degrees, printed r, theta, d, phi
ROWS = [
    ("Tetrahelix", "Tet", 6, 1, -120, .094, 131.810, -.516, 161.565),
    ("Tetratorus", "Tet", 3, 1, 0, .471, 70.529, 0, 90),
    ("Boxbeam", "Cub", 4, 1, -180, 0, 0, 1, 0),
    ("Staircase", "Cub", 4, 2, -180, 0, 0, .707, 0),
    ("Blockhelix", "Cub", 8, 2, -90, .236, 120, -.577, 144.736),
    ("Cubatorus", "Cub", 4, 2, 0, .500, 90, 0, 90),
    ("Octabeam", "Oct", 3, 5, -60, 0, 0, 1.155, 0),
    ("Octaspiky", "Oct", 6, 1, -120, .148, 146.443, -.603, 154.761),
    ("Octamedium", "Oct", 6, 2, -120, .163, 131.810, -.894, 161.565),
    ("Octagear", "Oct", 3, 1, 0, .408, 109.471, 0, 90),
    ("Treestar", "Oct", 3, 2, 0, .816, 70.529, 0, 90),
    ("Dodecabeam", "Dod", 5, 8, -108, 0, 0, 1.589, 0),
    ("Dodecadoubler", "Dod", 10, 1, -144, .113, 161.301, -.805, 164.550),
    ("The Alternater", "Dod", 10, 2, -144, .118, 149.520, -1.333, 170.306),
    ("Dodecashaft", "Dod", 10, 1, -72, .351, 129.657, -.543, 130.501),
    ("Dodecagear", "Dod", 5, 1, 0, .491, 116.565, 0, 90),
    ("Dodecacorkscrew", "Dod", 10, 2, -72, .546, 93.026, -1.095, 144.110),
    ("Dodecadonut", "Dod", 5, 2, 0, 1.286, 63.435, 0, 90),
    ("Pearlshaft", "Ico", 3, 13, -60, 0, 0, 1.589, 0),
    ("Quasi-planar", "Ico", 6, 8, 165, .049, 167.764, 1.294, 4.347),
    ("Two Strands", "Ico", 6, 1, 120, .137, 159.446, .499, 28.340),
    ("Slow Twist", "Ico", 6, 12, 120, .169, 124.309, 1.454, 11.641),
    ("Rock Candy", "Ico", 12, 2, 120, .204, 146.443, .830, 25.239),
    ("Icosa Tree Star", "Ico", 3, 1, 0, .304, 138.190, 0, 90),
    ("Icosacorkscrew", "Ico", 6, 8, -75, .512, 99.253, -1.037, 143.042),
    ("Planar point cluster", "Ico", 6, 2, 0, .562, 109.471, 0, 90),
    ("Big Icosacorkscrew", "Ico", 6, 8, 45, .803, 82.064, .756, 54.343),
    ("The Wheel", "Ico", 3, 12, 0, 2.080, 41.810, 0, 90),
]

def vertices(solid):
    """Canonical coordinates with exact unit edges, evaluated as floats."""
    g = (1 + sqrt(5)) / 2
    signs = list(product((-1, 1), repeat=3))
    if solid == "Tet":
        return [tuple(x / (2 * sqrt(2)) for x in p)
                for p in signs if p[0] * p[1] * p[2] == 1]
    if solid == "Cub":
        return [tuple(x / 2 for x in p) for p in signs]
    if solid == "Oct":
        return [tuple(s / sqrt(2) if i == j else 0.0 for i in range(3))
                for j in range(3) for s in (-1, 1)]
    if solid == "Dod":
        pts = list(signs)
        for s, t in product((-1, 1), repeat=2):
            pts += [(0, s / g, t * g), (s / g, t * g, 0),
                    (t * g, 0, s / g)]
        return [tuple(x * g / 2 for x in p) for p in pts]
    if solid == "Ico":
        pts = []
        for s, t in product((-1, 1), repeat=2):
            pts += [(0, s, t * g), (s, t * g, 0), (t * g, 0, s)]
        return [tuple(x / 2 for x in p) for p in pts]
    raise ValueError(solid)

def source_edge(solid):
    """Original THREE constructor edge length (Cube already has unit edges)."""
    w = sqrt(5)
    return {"Tet": sqrt(8/3), "Cub": 1.0, "Oct": sqrt(2),
            "Dod": 4/(sqrt(3)*(1+w)), "Ico": 4/sqrt(10+2*w)}[solid]

def inradius(solid):
    w = sqrt(5)
    return {"Tet": sqrt(6)/12, "Cub": .5, "Oct": sqrt(6)/6,
            "Dod": sqrt(250+110*w)/20,
            "Ico": sqrt(3)*(3+w)/12}[solid]

# Inserted below from the pinned repository's face numbering, after matching
# its Float32 vertices to these exact-coordinate vertex sets.
FACES = {'Tet': [[0, 3, 1], [3, 1, 2], [0, 3, 2], [0, 1, 2]],
 'Cub': [[7, 6, 5, 4],
         [2, 3, 0, 1],
         [7, 6, 2, 3],
         [5, 4, 0, 1],
         [7, 5, 3, 1],
         [6, 4, 2, 0]],
 'Oct': [[3, 5, 1],
         [5, 1, 2],
         [1, 2, 4],
         [3, 1, 4],
         [3, 4, 0],
         [2, 4, 0],
         [5, 2, 0],
         [3, 5, 0]],
 'Dod': [[17, 7, 3, 18, 12],
         [7, 18, 19, 13, 6],
         [13, 6, 4, 8, 14],
         [8, 14, 0, 10, 2],
         [0, 10, 9, 1, 16],
         [18, 12, 6, 14, 2],
         [3, 12, 10, 2, 16],
         [17, 3, 1, 16, 11],
         [4, 8, 0, 9, 15],
         [17, 7, 19, 11, 5],
         [19, 13, 4, 15, 5],
         [9, 1, 11, 15, 5]],
 'Ico': [[8, 9, 4],
         [9, 4, 10],
         [4, 10, 6],
         [4, 6, 2],
         [8, 4, 2],
         [9, 10, 11],
         [8, 9, 3],
         [8, 2, 1],
         [6, 2, 0],
         [10, 6, 5],
         [11, 3, 7],
         [3, 1, 7],
         [1, 0, 7],
         [0, 5, 7],
         [11, 5, 7],
         [9, 11, 3],
         [8, 3, 1],
         [2, 1, 0],
         [6, 0, 5],
         [10, 11, 5]]}
