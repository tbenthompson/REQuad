import numpy as np
import quadracheer as qc
from scipy.special import legendre
basis = legendre(16)
sing_pt = 0.5
gauss = {i:qc.gaussxw(i) for i in range(1, 17)}
gauss128 = qc.gaussxw(128)
func = lambda x: basis(x) * np.log(np.abs(x - sing_pt))
def integrate(x_min, x_max, q):
    x = q[0]
    w = q[1]
    new_w = w * ((x_max - x_min) / 2.0)
    new_x = x_min + (x_max - x_min) * ((x + 1.0) / 2.0)
    return sum(new_w * func(new_x))

def with_hole(h, q):
    return integrate(-1.0, sing_pt - h, q) + integrate(sing_pt + h, 1.0, q)

def donut(ho, hi, q):
    return integrate(sing_pt - ho, sing_pt - hi, q) + integrate(sing_pt + hi, sing_pt + ho, q)

n = 24
hs = 0.1 * (0.5 ** np.arange(0, n))

ests = []

for m in range(n):
    if m > 0:
        outer = hs[m - 1]
        previous = ests[0][-1]
    else:
        outer = 0.5
        previous = with_hole(0.5, gauss[16])
    addt = donut(outer, hs[m], gauss[max(3, 16 - 1 * m)])
    addt_ex = donut(outer, hs[m], gauss128)
    assert(abs(addt - addt_ex) < 1e-11)
    current = previous + addt
    ests.append([])
    ests[0].append(current)
    for i in range(1, m + 1):
        factor = 1.0 / (2.0 ** i - 1.0)
        prev_bad = ests[i-1][-2]
        prev_good = ests[i-1][-1]
        ests[i].append(factor * ((2.0 ** i) * prev_good - prev_bad))

exact = 0.0243191
print ests[-1][-1]
print ests[0][-1]
print abs(ests[-1][-1] - exact)
