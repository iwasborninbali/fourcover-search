"""Independent exact checks using only the Python standard library."""
from fractions import Fraction as F
import datetime
import hashlib
import json
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
INPUT = Path(sys.argv[1]).resolve() if len(sys.argv)>1 else ROOT / 'generator.json'
data = json.loads(INPUT.read_text())
k = 730033053750
Q = tuple(map(F, data['generator_Q']))
R = tuple(map(F, data['isogenous_curve']['generator']))

def add(P, Q, p, b):
    if P is None: return Q
    if Q is None: return P
    x, y = P
    u, v = Q
    if x == u and (y + v) % p == 0: return None
    slope = ((v-y)*pow(u-x, -1, p) if x != u else 3*x*x*pow(2*y, -1, p)) % p
    xx = (slope*slope-x-u) % p
    return xx, (slope*(x-xx)-y) % p

def mul(n, P, p, b):
    S = None
    while n:
        if n & 1: S = add(S, P, p, b)
        P = add(P, P, p, b)
        n //= 2
    return S

def points(p, b):
    roots = {}
    for y in range(p): roots.setdefault(y*y % p, []).append(y)
    return [None] + [(x,y) for x in range(p) for y in roots.get((x**3+b)%p, [])]

def reduce(P, p):
    x, y = P
    if x.denominator % p == 0:
        assert y.denominator % p == 0
        return None
    return x.numerator*pow(x.denominator,-1,p)%p, y.numerator*pow(y.denominator,-1,p)%p

def rational_add(P, Q):
    if P is None: return Q
    if Q is None: return P
    x,y=P
    u,v=Q
    if x==u and y==-v: return None
    s=(v-y)/(u-x) if x!=u else 3*x*x/(2*y)
    xx=s*s-x-u
    return xx,s*(x-xx)-y

def witness(P, b, ell, p):
    assert p > 3 and (432*b*b) % p != 0
    group = points(p,b)
    image = {mul(ell,T,p,b) for T in group}
    Pbar = reduce(P,p)
    assert Pbar not in image
    return {'ell':ell, 'p':p, 'order':len(group), 'image_size':len(image), 'Pbar':Pbar,
            'image_test':'enumerated every point and its ell-multiple'}

q_on = Q[1]**2-Q[0]**3 == k
r_on = R[1]**2-R[0]**3 == -k//27
assert q_on and r_on
# Standard degree-three map to b=-27k, followed by x/9,y/27.
x,y=Q
phi=((x**3+4*k)/(9*x*x), y*(x**3-8*k)/(27*x**3))
assert phi == R or (phi[0],-phi[1]) == R
xx,yy=R
dual=((xx**3+4*(-k//27))/(xx*xx), yy*(xx**3-8*(-k//27))/(xx**3))
assert dual[1]**2-dual[0]**3 == k
assert dual == rational_add(Q,rational_add(Q,Q))
# Torsion injects at good odd primes 17 and 31.
orders={str(p):len(points(p,k)) for p in (17,31)}
assert math.gcd(*orders.values()) == 1
assert {str(p):len(points(p,-k//27)) for p in (17,31)} == orders
witnesses={2:17,3:17,5:29,7:41,11:131,13:233,17:101,19:113,23:137,29:173,31:557,37:229,41:491,43:31,47:281}
certs=[witness(Q,k,ell,p) for ell,p in witnesses.items()]
# Once E=ZQ, an isogeny-image index divides 3. Exclude it directly.
isog_cert=None
for p in (17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101):
    if 432*(k//27)**2 % p:
        G=points(p,-k//27)
        if reduce(R,p) not in {mul(3,T,p,-k//27) for T in G}:
            isog_cert=witness(R,-k//27,3,p)
            break
assert isog_cert is not None
out={'checked_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
     'input_sha256':hashlib.sha256(INPUT.read_bytes()).hexdigest(),
     'Q_on_curve':q_on,'isog_R_on_curve':r_on,'phi_matches_sign':1 if phi==R else -1,
     'dual_phi_Q_equals_3Q':True,'isogenous_torsion_order':1,
     'good_prime_orders':orders,'torsion_order':1,'prime_nondisivibility':certs,
     'isog_3_nondivisibility':isog_cert,
     'displayed_height_index_ratio':167.3053526503390095/0.20377867626770474,
     'displayed_height_index_sqrt':math.sqrt(167.3053526503390095/0.20377867626770474)}
(ROOT/'independent-check.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
