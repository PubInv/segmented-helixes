#!/usr/bin/env python3
"""Reproduce the companion note's numerical checks; Python 3, no packages.

The catalog transcription is from MathematicJournal/StackingHelix.tex at
8be3839de2204b966c8aedb7c28f6a6f5c18076e. The baseline rotation is the
two-stage construction in js/segment_helix_math.js::adjoinPrism.
Checks are numerical; the note supplies the mathematical proof.
"""

from collections import Counter
from itertools import product
import math
import random


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def cross(a, b):
    return (a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0])


def scale(a, k):
    return tuple(k * x for x in a)


def minus(a, b):
    return tuple(x - y for x, y in zip(a, b))


def plus(a, b):
    return tuple(x + y for x, y in zip(a, b))


def unit(a):
    return scale(a, 1.0 / math.sqrt(dot(a, a)))


def clamp(x):
    return min(1.0, max(-1.0, x))


def transpose(a):
    return tuple(zip(*a))


def multiply(a, b):
    return tuple(tuple(dot(row, col) for col in transpose(b)) for row in a)


def apply(a, v):
    return tuple(dot(row, v) for row in a)


def matrix_error(a, b):
    return max(abs(x - y) for ra, rb in zip(a, b) for x, y in zip(ra, rb))


def rotation(axis, angle):
    """Rodrigues rotation; column-vector convention."""
    x, y, z = u = unit(axis)
    c, s = math.cos(angle), math.sin(angle)
    skew = ((0.0, -z, y), (z, 0.0, -x), (-y, x, 0.0))
    return tuple(tuple(c * (i == j) + (1.0 - c) * u[i] * u[j]
                       + s * skew[i][j] for j in range(3)) for i in range(3))


def perpendicular(n):
    k = min(range(3), key=lambda i: abs(n[i]))
    basis = tuple(float(i == k) for i in range(3))
    return unit(cross(n, basis))


def shortest(a, b):
    a, b = unit(a), unit(b)
    v = cross(a, b)
    sine = math.sqrt(dot(v, v))
    cosine = clamp(dot(a, b))
    if sine < 1e-14:
        return rotation(perpendicular(a), math.pi if cosine < 0 else 0.0)
    return rotation(scale(v, 1.0 / sine), math.atan2(sine, cosine))


def principal_angle(r):
    return math.acos(clamp((sum(r[i][i] for i in range(3)) - 1.0) / 2.0))


def half_angle(alpha, twist):
    return 2.0 * math.acos(clamp(math.cos(alpha / 2.0)
                               * abs(math.cos(twist / 2.0))))


def baseline(nb, nc, e):
    return multiply(shortest(e, nc), shortest(nb, scale(e, -1.0)))


def reference_phase(r0, a, nc):
    g = multiply(r0, transpose(a))
    v = perpendicular(nc)
    gv = apply(g, v)
    return math.atan2(dot(nc, cross(v, gv)), dot(v, gv))


def require_near(x, y, tol=2e-12):
    if abs(x - y) > tol:
        raise AssertionError(f"{x!r} != {y!r}; tolerance {tol}")


S2, S5, S6 = math.sqrt(2), math.sqrt(5), math.sqrt(6)
# name, solid, class, cos(alpha_j), twist degrees, exact cos(theta), printed theta
CASES = [
    ("Tetratorus", "Tet", "A", 1/3, 0, 1/3, 70.529),
    ("Tetrahelix", "Tet", "A", 1/3, -120, -2/3, 131.810),
    ("Cubatorus", "Cub", "A", 0, 0, 0, 90.000),
    ("Blockhelix", "Cub", "A", 0, -90, -1/2, 120.000),
    ("Octagear", "Oct", "A", -1/3, 0, -1/3, 109.471),
    ("Octaspiky", "Oct", "A", -1/3, -120, -5/6, 146.443),
    ("Treestar", "Oct", "S", 1/3, 0, 1/3, 70.529),
    ("Octamedium", "Oct", "S", 1/3, -120, -2/3, 131.810),
    ("Dodecagear", "Dod", "A", -1/S5, 0, -1/S5, 116.565),
    ("Dodecashaft", "Dod", "A", -1/S5, -72, (S5-15)/20, 129.657),
    ("Dodecadoubler", "Dod", "A", -1/S5, -144, -1/2-1/S5, 161.301),
    ("Dodecadonut", "Dod", "S", 1/S5, 0, 1/S5, 63.435),
    ("Dodecacorkscrew", "Dod", "S", 1/S5, -72, -1/2+1/S5, 93.026),
    ("The Alternater", "Dod", "S", 1/S5, -144, -3/4-S5/20, 149.520),
    ("Icosa Tree Star", "Ico", "A", -S5/3, 0, -S5/3, 138.190),
    ("Two Strands", "Ico", "A", -S5/3, 120, -(9+S5)/12, 159.446),
    ("Planar point cluster", "Ico", "S", -1/3, 0, -1/3, 109.471),
    ("Rock Candy", "Ico", "S", -1/3, 120, -5/6, 146.443),
    ("Big Icosacorkscrew", "Ico", "T", 1/3, 45, (S2-1)/3, 82.064),
    ("Icosacorkscrew", "Ico", "T", 1/3, -75, (S6-S2-2)/6, 99.253),
    ("Quasi-planar", "Ico", "T", 1/3, 165, -(2+S6+S2)/6, 167.764),
    ("The Wheel", "Ico", "F", S5/3, 0, S5/3, 41.810),
    ("Slow Twist", "Ico", "F", S5/3, 120, (S5-9)/12, 124.309),
]


def face_normals():
    phi = (1 + S5) / 2
    tet = [(1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1)]
    cub = [(s * (i == k), s * (j == k), s * (l == k))
           for k in range(3) for s in (-1, 1) for i, j, l in [(0, 1, 2)]]
    octa = list(product((-1, 1), repeat=3))
    dod = [(0, s, t * phi) for s, t in product((-1, 1), repeat=2)]
    dod += [(t * phi, 0, s) for s, t in product((-1, 1), repeat=2)]
    dod += [(s, t * phi, 0) for s, t in product((-1, 1), repeat=2)]
    ico = list(octa)
    for s, t in product((-1, 1), repeat=2):
        ico += [(0, s / phi, t * phi), (s / phi, t * phi, 0),
                (t * phi, 0, s / phi)]
    return {k: [unit(v) for v in vs] for k, vs in
            {"Tet": tet, "Cub": cub, "Oct": octa, "Dod": dod, "Ico": ico}.items()}


def verify_catalog():
    print("23 catalog checks: angles in degrees")
    print(f"{'Name':22s} {'Computed':>12s} {'Catalog':>10s} {'Abs. error':>12s}")
    worst = 0.0
    for name, solid, cls, ca, twist, exact, printed in CASES:
        value = (1 + ca) * (1 + math.cos(math.radians(twist))) / 2 - 1
        require_near(value, exact)
        theta = math.degrees(half_angle(math.acos(ca), math.radians(twist)))
        require_near(theta, math.degrees(math.acos(clamp(exact))), tol=1e-10)
        error = abs(theta - printed)
        if error >= 0.0005:
            raise AssertionError(f"{name}: failed printed-precision comparison")
        worst = max(worst, error)
        print(f"{name:22s} {theta:12.7f} {printed:10.3f} {error:12.7f}")
    print(f"PASS: 23/23; largest error {worst:.7f} degrees (< 0.0005).")


def verify_geometry():
    expected = {
        "Tet": {-1/3: 3},
        "Cub": {0: 4, -1: 1},
        "Oct": {1/3: 3, -1/3: 3, -1: 1},
        "Dod": {1/S5: 5, -1/S5: 5, -1: 1},
        "Ico": {S5/3: 3, 1/3: 6, -1/3: 6, -S5/3: 3, -1: 1},
    }
    comparisons = 0
    for solid, ns in face_normals().items():
        counts = Counter(round(dot(ns[0], n), 9) for n in ns[1:])
        target = Counter({round(s, 9): count for s, count in expected[solid].items()})
        if counts != target:
            raise AssertionError((solid, counts, target))
        for nb in ns:
            for nc in ns:
                s = dot(nb, nc)
                if abs(s) > 1 - 1e-12:
                    continue
                e = unit(minus(nc, nb))
                r0 = baseline(nb, nc, e)
                a = shortest(nb, scale(nc, -1))
                if matrix_error(r0, a) > 3e-12:
                    raise AssertionError("Face-center baseline was not shortest alignment")
                comparisons += 1
        print(f"PASS: {solid} face-normal multiplicities.")
    print(f"PASS: {comparisons} ordered nonparallel face-center baselines have zero phase.")


def verify_rotations():
    rng = random.Random(20261001)
    for _ in range(1000):
        alpha = rng.uniform(0.001, math.pi - 0.001)
        twist = rng.uniform(-4 * math.pi, 4 * math.pi)
        r = multiply(rotation((0, 0, 1), twist), rotation((1, 0, 0), alpha))
        require_near(principal_angle(r), half_angle(alpha, twist), tol=3e-10)
    print("PASS: 1000 perpendicular-axis matrix/half-angle comparisons.")

    for _ in range(250):
        nb = unit(tuple(rng.uniform(-1, 1) for _ in range(3)))
        nc = unit(tuple(rng.uniform(-1, 1) for _ in range(3)))
        e = (0, 0, 1)
        a = shortest(nb, scale(nc, -1))
        r0 = baseline(nb, nc, e)
        gamma = reference_phase(r0, a, nc)
        if matrix_error(r0, multiply(rotation(nc, gamma), a)) > 3e-12:
            raise AssertionError("Reference-phase factorization failed")
        twist = rng.uniform(-math.pi, math.pi)
        r = multiply(rotation(nc, twist), r0)
        alpha = math.acos(clamp(-dot(nb, nc)))
        require_near(principal_angle(r), half_angle(alpha, twist + gamma), tol=3e-10)
    print("PASS: 250 general two-stage baselines with computed reference phases.")

    nb, nc, e = unit((1, 0, -1)), unit((0, 1, 1)), (0, 0, 1)
    r0 = baseline(nb, nc, e)
    a = shortest(nb, scale(nc, -1))
    gamma = reference_phase(r0, a, nc)
    require_near(gamma, math.acos(2*S2/3))
    require_near(principal_angle(r0), math.acos((2*S2-1)/4))
    print(f"PASS: noncoplanar example: alpha = 60 degrees, phase = "
          f"{math.degrees(gamma):.7f}, full rotation = "
          f"{math.degrees(principal_angle(r0)):.7f} degrees.")


def verify_straight_joints():
    cases = [
        ("Boxbeam", (0, 0, -1), (0, 0, 1), -180, 180),
        ("Staircase", (1, 0, 0), (0, 1, 0), -180, 180),
        ("Octabeam", (0, 0, -1), (0, 0, 1), -60, 60),
        ("Dodecabeam", (0, 0, -1), (0, 0, 1), -108, 108),
        ("Pearlshaft", (0, 0, -1), (0, 0, 1), -60, 60),
    ]
    for name, nb, nc, twist, expected in cases:
        e = unit(minus(nc, nb))
        r = multiply(rotation(nc, math.radians(twist)), baseline(nb, nc, e))
        require_near(math.degrees(principal_angle(r)), expected, tol=2e-6)
        # Unit inradius; T sends B=nb to C=nc. Check T(C)-C == C-B.
        b, c = nb, nc
        t = minus(c, apply(r, b))
        d = plus(apply(r, c), t)
        if math.sqrt(dot(minus(minus(d, c), minus(c, b)),
                         minus(minus(d, c), minus(c, b)))) > 3e-12:
            raise AssertionError(f"{name}: joint points were not collinear/equally spaced")
        print(f"PASS: {name}: full rotation {expected} degrees; straight joint orbit.")


def main():
    verify_catalog()
    print()
    verify_geometry()
    verify_rotations()
    verify_straight_joints()
    print("\nAll checks passed. Numerical checks supplement the proof; they do not replace it.")


if __name__ == "__main__":
    main()
