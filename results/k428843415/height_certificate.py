"""Exact rational index bound: Voutier--Yabuta 1.2(a) and 1.3.

Fresh implementation. ln(n) is enclosed by a dyadic enclosure of the atanh
argument, its positive power series, and an explicit geometric tail bound.
No floating-point operation participates in any proof comparison.
"""
from fractions import Fraction as F
from pathlib import Path
import json

ROOT = Path(__file__).resolve().parent
data = json.loads((ROOT/'generator.json').read_text())
B = 428843415
x = F(data['generator_Q'][0])

def log_bounds(n):
    assert isinstance(n,int) and n>0
    e = n.bit_length()-1
    def series(t):
        # t is in [0,1/3], except for a tiny outward dyadic rounding.
        assert 0 <= t < F(1,2)
        term=t
        s=F(0)
        for j in range(96):
            s += 2*term/F(2*j+1)
            term *= t*t
        tail = 2*term/(193*(1-t*t))
        return s,s+tail
    def enclosure(t):
        grid=2**128
        a=(t*grid).__floor__()
        return series(F(a,grid))[0], series(F(a+1,grid))[1]
    l2,u2=enclosure(F(1,3))
    # n=2^e*m, and (m-1)/(m+1)=(n-2^e)/(n+2^e).
    l,u=enclosure(F(n-2**e,n+2**e))
    l,u=l+e*l2,u+e*u2
    scale=10**30
    return F((l*scale).__floor__(),scale),F((u*scale).__ceil__(),scale)

L2,U2=log_bounds(2)
L3,U3=log_bounds(3)
LB,UB=log_bounds(B)
A=max(abs(x.numerator),x.denominator)
LA,UA=log_bounds(A)
# h_PAR(P)=2*h_VY(P); Theorem 1.2(a), positive b.
low=LB/3-F(4,3)*U2-F(3,2)*U3-F(12,1000)
low_hi=UB/3-F(4,3)*L2-F(3,2)*L3-F(12,1000)
# Uniform lower bound in Theorem 1.3 for (h_naive/2 - h_VY).
up=UA+UB/3+F(598,1000)
up_lo=LA+LB/3+F(598,1000)
assert low > 0
bound=0
while (bound+1)**2*low < up:
    bound+=1
assert up < (bound+1)**2*low
assert bound == 7
out=dict(theorems=['1.2(a)','1.3'], normalization='H_PAR=2*h_VY',
         A=str(A), log_B=[str(LB),str(UB)], log_A=[str(LA),str(UA)],
         H_all_nontorsion_lower=[str(low),str(low_hi)],
         H_Q_upper=[str(up_lo),str(up)],
         exact_upper_over_lower=str(up/low),
         exact_upper_less_than_64_lower=up<64*low,
         index_at_most=bound, primes_to_exclude=[2,3,5,7],
         decimal_orientation_only=dict(lower=float(low),upper=float(up),ratio=float(up/low)),
         method='96 atanh terms, dyadic argument enclosure 2^-128, exact positive tail bound, outward 10^-30 rounding')
(ROOT/'height-certificate.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
