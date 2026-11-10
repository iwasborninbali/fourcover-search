"""Independent standard-library arithmetic for the 428843415 audit."""
from fractions import Fraction as F
from pathlib import Path
import datetime
import hashlib
import json
import math

ROOT = Path(__file__).resolve().parent
INPUT = ROOT / 'generator.json'
data = json.loads(INPUT.read_text())
B = 428843415
Q = tuple(F(s) for s in data['generator_Q'])
R = tuple(F(s) for s in data['isogenous_curve']['generator'])

def plus(P, T, modulus=None):
    if P is None:
        return T
    if T is None:
        return P
    x, y = P
    u, v = T
    mod = (lambda z: z % modulus) if modulus else (lambda z: z)
    if x == u and mod(y + v) == 0:
        return None
    top, bottom = (v-y, u-x) if x != u else (3*x*x, 2*y)
    s = mod(top * pow(bottom, -1, modulus)) if modulus else top / bottom
    a = mod(s*s-x-u)
    return a, mod(s*(x-a)-y)

def times(n, P, modulus=None):
    ans = None
    while n:
        if n % 2:
            ans = plus(ans, P, modulus)
        P = plus(P, P, modulus)
        n //= 2
    return ans

def group(b, p):
    squares = {}
    for y in range(p):
        squares.setdefault(y*y % p, []).append(y)
    return [None] + [(x, y) for x in range(p) for y in squares.get((x**3+b) % p, [])]

def reduction(P, p):
    if P[0].denominator % p == 0:
        assert P[1].denominator % p == 0
        return None
    return tuple(a.numerator * pow(a.denominator, -1, p) % p for a in P)

def certificate(P, b, ell, p):
    assert p > 3 and (432*b*b) % p
    G = group(b, p)
    image = {times(ell, T, p) for T in G}
    bar = reduction(P, p)
    assert bar not in image, (b, ell, p, bar)
    return dict(ell=ell, p=p, order=len(G), image_size=len(image), reduction=bar,
                outside_multiplication_image=True)

assert B == 3*5*7*11*13**5
assert Q[1]**2 == Q[0]**3 + B
assert R[1]**2 == R[0]**3 - 27*B
x, y = Q
phi = ((x**3+4*B)/x**2, y*(x**3-8*B)/x**3)
assert phi == R
u, v = R
dual = ((u**3-108*B)/(9*u**2), v*(u**3+216*B)/(27*u**3))
assert dual == times(3, Q)
assert dual[1]**2 == dual[0]**3+B

# An identity of polynomials verifies the map on the entire affine domain x != 0.
# (X+B)(X-8B)^2 - (X+4B)^3 + 27B X^2 = 0, X=x^3.
def polymul(a, b):
    ans = [0]*(len(a)+len(b)-1)
    for i, v in enumerate(a):
        for j, w in enumerate(b):
            ans[i+j] += v*w
    return ans
a = polymul([B, 1], polymul([-8*B, 1], [-8*B, 1]))
b = polymul([4*B, 1], polymul([4*B, 1], [4*B, 1]))
assert [a[i]-b[i]+(27*B if i == 2 else 0) for i in range(4)] == [0]*4

# Find a coprime pair of actual good-reduction group orders, without assuming
# the torsion routine or reusing the sister curve's choice of primes.
small_primes = [17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101]
orders = {p:len(group(B,p)) for p in small_primes if (432*B*B)%p}
torsion_pair = next((p,q) for p in orders for q in orders if p<q and math.gcd(orders[p], orders[q]) == 1)
torsion_orders = {p:orders[p] for p in torsion_pair}
partner_orders = {p:len(group(-27*B,p)) for p in torsion_pair}
assert partner_orders == torsion_orders

witness_primes = {2:17,3:17,5:29,7:41,11:131,13:233,17:101,19:113,
                  23:137,29:173,31:43,37:439,41:983,43:31,47:281}
certs = [certificate(Q, B, ell, p) for ell,p in witness_primes.items()]
partner_cert = next(certificate(R,-27*B,3,p) for p in small_primes
                    if (432*B*B)%p and reduction(R,p) not in {times(3,T,p) for T in group(-27*B,p)})
out = dict(checked_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
           input=str(INPUT), input_sha256=hashlib.sha256(INPUT.read_bytes()).hexdigest(),
           Q_on_curve=True, R_on_partner=True, phi_exactly_matches_stored_R=True,
           isogeny_polynomial_identity=True, dual_R_equals_3Q=True,
           good_reduction_orders=torsion_orders, partner_orders=partner_orders,
           prime_nondivisibility=certs, partner_3_nondivisibility=partner_cert)
(ROOT/'exact-check.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
