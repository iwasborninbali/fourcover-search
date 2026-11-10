# Saturation of the found point on y^2 = x^3 + 730033053750 (Sage: Cremona-Prickett-Siksek height bounds).
import json, time
from sage.all import EllipticCurve, QQ
r = next(json.loads(l) for l in open("verify-M.log") if l.startswith("{") and json.loads(l).get("verified"))
v = r["verified"][0]
E = EllipticCurve(QQ, [0, 0, 0, 0, 730033053750]); Q = E([QQ(v["Q"][0]), QQ(v["Q"][1])])
print("Q on E, order:", Q.order(), "torsion:", E.torsion_order(), flush=True)
t0 = time.time()
lam = E.height_function().min(0.1, 20) if hasattr(E, "height_function") else None
print("CPS-type lower bound for heights of nontorsion points (Sage height_function().min):", lam, "seconds", round(time.time() - t0, 1), flush=True)
t0 = time.time()
sat = E.saturation([Q], verbose=True)
print("saturation result (points, index, regulator):", sat, "seconds", round(time.time() - t0, 1), flush=True)
