# Relative four-descent above one 2-cover v^2 = R(t) of E, reduction of each soluble four-cover, and exact verification.
# Four-descent: C. D. Long's elliptic-descent (https://github.com/long-algorithms/elliptic-descent, MIT; tested at
# commit 882a0dac), proof=True. Set ELLIPTIC_DESCENT_SRC to the src/ directory of that checkout.
# Reduction (heuristic, self-verifying): U in c^{-1} with (xi) = b c^2, LLL against
# T_xi(U) = sum_sigma |sigma(xi)| |sigma(U)|^2. Every point reported downstream is verified exactly here.
# Usage (run with Sage's Python):
#   WORLD_JSON=world.json sage -python search4.py COVER 0              -> one JSON line per soluble four-cover, with
#                                                                        the reduced quadrics Q2, Q3 for elkies4q
#   WORLD_JSON=world.json sage -python search4.py COVER 0 points.txt   -> exact check of points "label u0 u1 u2 u3 1",
#                                                                        mapped to E (the 'verified' field)
#   WORLD_JSON=world.json sage -python search4.py COVER B              -> naive box search |u0|,|u1|,|u2| <= B (small B)
import json, sys, os, time
from math import isqrt
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.environ.get('ELLIPTIC_DESCENT_SRC', os.path.join(HERE, 'src')))
from sage.all import QQ, ZZ, PolynomialRing, NumberField, Matrix, RealField, prod
from elliptic_descent.relative import relative_four_descent

# PARI ell2cover basis with exact maps (ascending coefficients), written by export_world.gp
W3 = json.load(open(os.environ['WORLD_JSON']))
COVERS = {cv['label']: tuple(cv['R']) for cv in W3['covers']}
which, B = sys.argv[1], int(sys.argv[2])
F = COVERS[which]
log = {'cover': which, 'ascending': F, 'box': B, 'candidates': []}
t0 = time.time()
R = relative_four_descent(F, proof=True)
log['descent'] = {'seconds': time.time() - t0, 'classes': R.total_model_classes,
                  'proof': str(R.arithmetic_proof), 'assumptions': [str(a) for a in R.assumptions]}

Pt = PolynomialRing(QQ, 't'); t = Pt.gen()
f = sum(F[i] * t**i for i in range(5))
K = NumberField(f, 'a'); a = K.gen()
nbits = lambda q: QQ(q).numerator().nbits() + QQ(q).denominator().nbits()
S = PolynomialRing(QQ, ['u0', 'u1', 'u2', 'u3']); us = S.gens()
ST = PolynomialRing(S, 'T'); T = ST.gen()
fT = ST(list(F))

# Optional third argument: a file of points "label u0 u1 u2 num den" (u3 = num/den) to verify exactly.
VERIFY = None
if len(sys.argv) > 3:
    VERIFY = [tuple(int(z) for z in line.split()) for line in open(sys.argv[3]) if line.strip()]
COV = next(cv for cv in W3['covers'] if cv['label'] == which)
from sage.all import EllipticCurve
E = EllipticCurve(QQ, W3['ainvs']); P0 = E([QQ(z) for z in W3['known_point']]) if W3.get('known_point') else None

def verify_point(label, qs, Q2, Q3, p):
    """Point on the four-cover -> t on the two-cover v^2 = R(t) -> frozen map to E; every step exact."""
    _, u0, u1, u2, num, den = p
    u = [QQ(u0), QQ(u1), QQ(u2), QQ(num) / QQ(den)]
    out = {'u': [str(z) for z in u], 'on_four_cover': bool(Q2(u) == 0 and Q3(u) == 0)}
    q0, q1 = QQ(qs[0](u)), QQ(qs[1](u))
    if q1 == 0:
        out['t'] = 'infinity'; return out
    tt = -q0 / q1; Rt = sum(QQ(F[i]) * tt**i for i in range(5)); out['t'] = str(tt)
    out['R_t_is_square'] = bool(Rt.is_square())
    if not out['R_t_is_square'] or Rt == 0:
        return out
    v = Rt.sqrt()
    Nt = sum(QQ(c) * tt**i for i, c in enumerate(COV['N'])); Tt = sum(QQ(c) * tt**i for i, c in enumerate(COV['T']))
    a1, a3 = QQ(W3['ainvs'][0]), QQ(W3['ainvs'][2])
    x = Nt / v**2; y = (-(a1 * x + a3) + Tt / v**3) / 2
    out['on_E'] = bool(E.is_on_curve(x, y))
    if out['on_E']:
        Q = E([x, y]); out['Q'] = [str(x), str(y)]
        out['height_Q'] = float(Q.height())
        if P0 is not None: out['regulator_P_Q'] = float(E.regulator_of_points([P0, Q]))
    return out

CALIB = [int(k) for k in os.environ.get('CALIB', '').split(',') if k]

def calibrate(k, xi, newb):
    """Exact preimage of kP: t with pi(t) = kP on v^2 = R(t), then U with xi*U^2 = mu*(theta - t);
    returns the primitive integer coordinates of U in the reduced basis (None if kP is not above this cover)."""
    Q = k * P0; xq, yq = Q[0], Q[1]
    Nt = sum(QQ(c) * t**i for i, c in enumerate(COV['N'])); Tt = sum(QQ(c) * t**i for i, c in enumerate(COV['T']))
    out = {'k': k, 'height': float(Q.height()), 'u': None}
    for tt, _ in (Nt - xq * f).roots(QQ):
        Rt = f(tt)
        if not Rt.is_square() or Rt == 0: continue
        gam = (a - tt) / xi
        ps = sorted(set([2] + [p for p, _ in QQ(gam.norm()).factor()] + [p for p, _ in QQ(f.leading_coefficient()).factor()]))
        for mask in range(2 ** (len(ps) + 1)):
            mu = prod([ps[i] for i in range(len(ps)) if mask >> i & 1], QQ(1)) * (-1 if mask >> len(ps) & 1 else 1)
            if (mu * gam).is_square():
                U = (mu * gam).sqrt()
                M = Matrix(QQ, [b.list() for b in newb])
                u = M.solve_left(Matrix(QQ, [U.list()]))[0]
                den = ZZ(1)
                for z in u: den = den.lcm(QQ(z).denominator())
                ui = [ZZ(z * den) for z in u]; g = ZZ(0)
                for z in ui: g = g.gcd(z)
                ui = [z // g for z in ui]
                out.update({'t': str(tt), 'u': [str(z) for z in ui], 'log_max_u': float(RealField(53)(max(abs(z) for z in ui)).log())})
                return out
    return out

def search(Q2, Q3, B):
    """All primitive (u0,u1,u2,u3), u0..u2 in box, u3 rational root of Q2, with Q3 = 0."""
    mons = lambda q: {tuple(int(e) for e in m.exponents()[0]): int(c) for c, m in zip(q.coefficients(), q.monomials())}
    m2, m3 = mons(Q2), mons(Q3)
    A = m2.get((0, 0, 0, 2), 0)
    lin = [m2.get(tuple(1 if k in (i, 3) else 0 for k in range(4)), 0) for i in range(3)]
    quad = {k: v for k, v in m2.items() if k[3] == 0}
    q3 = list(m3.items())
    found = []
    for x0 in range(0, B + 1):
        for x1 in range(-B, B + 1):
            for x2 in range(-B, B + 1):
                if x0 == 0 and (x1 < 0 or (x1 == 0 and x2 <= 0)):
                    continue
                x = (x0, x1, x2)
                Bc = lin[0] * x0 + lin[1] * x1 + lin[2] * x2
                C = 0
                for (e0, e1, e2, _), c in quad.items():
                    C += c * x0**e0 * x1**e1 * x2**e2
                roots = []
                if A == 0:
                    if Bc != 0: roots = [Fr(-C, Bc)]
                else:
                    D = Bc * Bc - 4 * A * C
                    if D >= 0:
                        s = isqrt(D)
                        if s * s == D:
                            roots = {Fr(-Bc + s, 2 * A), Fr(-Bc - s, 2 * A)}
                for r in roots:
                    v = (x0, x1, x2, r)
                    if sum(c * v[0]**e[0] * v[1]**e[1] * v[2]**e[2] * v[3]**e[3] for e, c in q3) == 0:
                        found.append(v)
    return found

for c in R.candidates:
    st = str(c.local.status)
    rec = {'label': c.label, 'local_status': st}
    log['candidates'].append(rec)
    if 'INSOLUBLE' in st:
        continue
    xi = K([QQ(str(z)) for z in c.model.xi])
    fac = K.ideal(xi).factor()
    cc = prod([P**(e // 2) for P, e in fac], K.ideal(1))
    basis = (cc**-1).basis()
    # Long's xi can carry a huge unit part (hundreds to thousands of digits). The embeddings then
    # span ~2*xb bits and evaluating them cancels ~xb bits, so the precision and the integer scale
    # of the Gram matrix must grow with xb (a fixed 600-bit / 2^120 setting left the model unreduced).
    xb = max(nbits(z) for z in xi.list()); bb = max(nbits(z) for b in basis for z in b.list())
    scale = 2 * xb + 4 * bb + 1000; prec = scale + 2 * xb + 4 * bb + 1000
    places = K.places(prec=prec)
    G = [[0] * 4 for _ in range(4)]; lognorm = 0
    for s in places:
        w = abs(s(xi)); mult = 1 if 'Real' in str(s.codomain()) else 2
        lognorm += mult * w.log()
        v = [s(b) for b in basis]
        for i in range(4):
            for j in range(4):
                G[i][j] += mult * w * (v[i] * v[j].conjugate()).real()
    err = abs(lognorm - RealField(prec)(abs(xi.norm())).log())
    assert err < 2.0**-100, f'precision check failed: |log N(xi)| error {err}'
    mx = max(abs(G[i][i]) for i in range(4))
    Gz = Matrix(ZZ, 4, 4, [[(G[i][j] / mx * 2**scale).round() for j in range(4)] for i in range(4)])
    rec['xi_bits'] = int(xb); rec['prec_bits'] = int(prec)
    Uni = Gz.LLL_gram()
    newb = [sum(Uni[i, j] * basis[i] for i in range(4)) for j in range(4)]
    coords = [b.list() for b in newb]  # power-basis coordinates
    Uelt = sum(sum(coords[j][k] * us[j] for j in range(4)) * T**k for k in range(4))
    prodT = (ST(xi.list()) * Uelt**2) % fT
    qs = [prodT[k] for k in range(4)]
    Q2, Q3 = S(qs[2]), S(qs[3])
    l2 = ZZ(1)
    for x in Q2.coefficients(): l2 = l2.lcm(QQ(x).denominator())
    l3 = ZZ(1)
    for x in Q3.coefficients(): l3 = l3.lcm(QQ(x).denominator())
    Q2, Q3 = Q2 * l2, Q3 * l3
    rec['reduced_coeff_digits'] = max(len(str(abs(ZZ(x)))) for x in list(Q2.coefficients()) + list(Q3.coefficients()))
    # Reduced quadrics for the compiled search (elkies4q.c): coefficient of u_i u_j for i <= j.
    qc = lambda q: [int(q.monomial_coefficient(us[i] * us[j])) for i in range(4) for j in range(i, 4)]
    rec['Q2'], rec['Q3'] = qc(Q2), qc(Q3)
    if CALIB:  # pull back kP to this four-cover (when it lies above it) and report its reduced coordinates
        rec['calibration'] = [calibrate(k, xi, newb) for k in CALIB]
        print(json.dumps(rec), flush=True)
        continue
    if VERIFY is not None:  # exact verification of points reported by elkies4q / elkies4 (no box search)
        rec['verified'] = [verify_point(c.label, qs, Q2, Q3, p) for p in VERIFY if p[0] == int(c.label)]
        print(json.dumps(rec), flush=True)
        continue
    t1 = time.time()
    pts = search(Q2, Q3, B)
    rec['search_seconds'] = time.time() - t1
    rec['points'] = []
    for v in pts[:20]:
        sub = dict(zip(us, [QQ(str(z)) for z in v]))
        q0, q1 = QQ(qs[0].subs(sub)), QQ(qs[1].subs(sub))
        Ue = sum(QQ(str(v[j])) * newb[j] for j in range(4))
        rec['points'].append({'u': [str(z) for z in v], 'q0': str(q0), 'q1': str(q1),
                              'normU': str(Ue.norm()), 'xi_norm': str(xi.norm())})
    print(json.dumps(rec), flush=True)

json.dump(log, open(os.path.join(os.environ.get('OUT_DIR', '.'), f"search4-{W3['world']}-{which}-B{B}.json"), 'w'), indent=1)
print('done', time.time() - t0, flush=True)
