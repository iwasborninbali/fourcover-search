# Image of the generator under the 3-isogeny to y^2 = x^3 - 11578772205, and saturation there.
import json
from sage.all import EllipticCurve, QQ
r = next(json.loads(l) for l in open("verify-M.log") if l.startswith("{") and json.loads(l).get("verified"))
v = r["verified"][0]
E = EllipticCurve(QQ, [0, 0, 0, 0, 428843415]); Q = E([QQ(v["Q"][0]), QQ(v["Q"][1])])
print("k factorization:", QQ(428843415).factor())
for phi in E.isogenies_prime_degree(3):
    E2 = phi.codomain(); Em = E2.minimal_model(); iso = E2.isomorphism_to(Em)
    R = iso(phi(Q))
    print("codomain minimal model:", Em.a_invariants(), "conductor", Em.conductor())
    sat = Em.saturation([R])
    print("image point:", R.xy())
    print("height", R.height(), "saturation index", sat[1], "regulator", sat[2])
    if sat[1] != 1: print("saturated generator:", sat[0][0].xy(), "height", sat[0][0].height())
