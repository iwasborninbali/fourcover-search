# Exact solver for the DENSE slices printed by elkies4.c (run with Sage's Python).
# A slice {x0 v0 + x1 v1 + w} spans a projective plane; it meets the four-cover D: Q2 = Q3 = 0 in at most four points.
# We find every rational point of D on that plane: z = a v0 + b v1 + c w with Q2(z) = Q3(z) = 0, (a:b:c) in P^2(Q).
# Usage: python dense_solve.py coeffs.txt elkies-output.txt   -> prints POINT lines (primitive, exactly checked)
import sys
from sage.all import QQ, ZZ, PolynomialRing, gcd

c = [int(z) for z in open(sys.argv[1]).read().split()]
q2, q3 = c[:10], c[10:]
R = PolynomialRing(QQ, ['a', 'b', 'cc']); a, b, cc = R.gens()

def Q(coef, z):
    k = 0; r = 0
    for i in range(4):
        for j in range(i, 4):
            r += coef[k] * z[i] * z[j]; k += 1
    return r

def primitive(v):
    v = [QQ(x) for x in v]
    d = ZZ(1)
    for x in v: d = d.lcm(x.denominator())
    w = [ZZ(x * d) for x in v]; g = ZZ(0)
    for x in w: g = g.gcd(x)
    w = [x // g for x in w]
    s = next((x for x in w if x), 1)
    return tuple(-x if s < 0 else x for x in w)

found, slices = set(), set()
for line in open(sys.argv[2]):
    if not line.startswith('DENSE'):
        continue
    n = [int(z) for z in line.split()[1:]]
    key = tuple(n)
    if key in slices:
        continue
    slices.add(key)
    v0, v1, w = n[0:4], n[4:8], n[8:12]
    # exact basis of the span of v0, v1, w (w = 0 near a rational point: the slice is a line through it)
    from sage.all import matrix
    basis = [list(r) for r in matrix(ZZ, [v0, v1, w]).row_module().basis()]
    sols = []
    if len(basis) == 1:
        sols.append(basis[0])
    elif len(basis) == 2:   # a line in P^3: common roots of two binary quadratic forms
        S2 = PolynomialRing(QQ, ['s', 'r']); sv, rv = S2.gens()
        z = [sv * basis[0][i] + rv * basis[1][i] for i in range(4)]
        g = gcd(Q(q2, z), Q(q3, z))
        for fac, _ in (g.factor() if g.degree() > 0 else []):
            if fac.degree() == 1:
                al, ga = fac.monomial_coefficient(sv), fac.monomial_coefficient(rv)
                sols.append([ga * basis[0][i] - al * basis[1][i] for i in range(4)])
    else:                   # a plane in P^3: it meets D in at most four points; no common factor is possible
        e0, e1, e2 = basis
        z = [a * e0[i] + b * e1[i] + cc * e2[i] for i in range(4)]
        f2, f3 = Q(q2, z), Q(q3, z)
        if f2.subs(a=0, cc=0) == 0 and f3.subs(a=0, cc=0) == 0:
            sols.append(e1)
        res = f2.resultant(f3, b)
        if res == 0:
            print('DEGENERATE', line.strip()); continue
        Sb = PolynomialRing(QQ, 'bb')
        for fac, _ in res.factor():
            if fac.degree() != 1:
                continue
            alpha, gamma = fac.monomial_coefficient(a), fac.monomial_coefficient(cc)
            av, cv = gamma, -alpha
            g = gcd(f2.subs(a=av, cc=cv).univariate_polynomial(Sb), f3.subs(a=av, cc=cv).univariate_polynomial(Sb))
            for r, _ in (g.roots(QQ) if g.degree() > 0 else []):
                sols.append([av * e0[i] + r * e1[i] + cv * e2[i] for i in range(4)])
    for s in sols:
        p = primitive(s)
        if any(p) and Q(q2, p) == 0 and Q(q3, p) == 0 and p not in found:
            found.add(p); print('POINT', *p, flush=True)
print(f'DENSE_DONE slices={len(slices)} points={len(found)}')
