"""A source-theorem based index certificate, using exact rational log intervals.

Voutier--Yabuta, arXiv:1305.6560v4, Theorems 1.2(c), 1.3.
Their canonical height is HALF the PARI/Sage/eclib convention here.
Local Tamagawa hypotheses are checked in certified_rank.out.
"""
from fractions import Fraction as F
import json
import sys
from pathlib import Path

ROOT=Path(__file__).resolve().parent
k=730033053750
INPUT=Path(sys.argv[1]).resolve() if len(sys.argv)>1 else ROOT/'generator.json'
data=json.loads(INPUT.read_text())
x=F(data['generator_Q'][0])

def atanh_log(z,N=100):
    assert 0<=z<1
    s=2*sum((z**(2*j+1)/F(2*j+1) for j in range(N)),F(0))
    err=2*z**(2*N+1)/(F(2*N+1)*(1-z*z))
    return s,s+err

L2,U2=atanh_log(F(1,3))
def log_interval(n):
    assert isinstance(n,int) and n>=1
    shift=n.bit_length()-1
    t=F(n,2**shift)
    lo,hi=atanh_log((t-1)/(t+1))
    lo,hi=lo+shift*L2,hi+shift*U2
    scale=10**50
    return F((lo*scale).__floor__(),scale),F((hi*scale).__ceil__(),scale)

Lk,Uk=log_interval(k)
L3,U3=log_interval(3)
La,Ua=log_interval(max(abs(x.numerator),x.denominator))
# h_PAR(P) > (log k)/9 -(2/3)log 2 -(1/3)log 3 -0.008.
lower=Lk/9-F(2,3)*U2-U3/3-F(8,1000)
lower_upper=Uk/9-F(2,3)*L2-L3/3-F(8,1000)
# Theorem 1.3: h_PAR(Q) < h_naive(Q)+(log k)/3+0.598.
upper=Ua+Uk/3+F(598,1000)
upper_lower=La+Lk/3+F(598,1000)
assert lower>0
assert lower>F(21988,10000) and upper<F(173433,1000)
assert upper<81*lower
assert upper/lower<81
assert upper/lower>64
primes=[2,3,5,7,11,13]
factor={}
n=k
for p in primes:
    e=0
    while n%p==0: n//=p;e+=1
    factor[p]=e
assert n==1 and all(e<6 for e in factor.values())
assert factor[5]==4 and (k//5**4)%5==1
assert factor[7]==5 and factor[11]==factor[13]==1
out={'theorem_source':'https://arxiv.org/abs/1305.6560v4',
     'theorems':['1.2(c)','1.3'],
     'factorization':factor,
     'height_convention':'PARI height = 2 * Voutier-Yabuta canonical height',
     'tamagawa_hypothesis':'For p>3: c5=3; c7=c11=c13=1; good elsewhere. Confirmed by GP minimal local reduction.',
     'all_nontorsion_PAR_height_lower_interval':[str(lower),str(lower_upper)],
     'Q_PAR_height_upper_interval':[str(upper_lower),str(upper)],
     'lower_decimal_for_orientation':float(lower),'upper_decimal_for_orientation':float(upper),
     'index_squared_upper_decimal_for_orientation':float(upper/lower),
     'exact_comparison_upper_less_than_81_lower':upper<81*lower,
     'index_at_most':8,'primes_requiring_nondisivibility':[2,3,5,7],
     'log_interval_method':'ln(n)=shift*ln(2)+2*sum z^(2j+1)/(2j+1), j=0..99; positive tail bounded by 2z^201/(201*(1-z^2)); rational arithmetic throughout'}
(ROOT/'height-certificate.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:v for k,v in out.items() if not k.endswith('_interval')},indent=2))
