import numpy as np
import quadracheer as qc
from scipy.special import legendre
basis = legendre(16)
sing_pt = 0.5
gauss = {i:qc.gaussxw(i) for i in range(1, 512, 8)}
gauss_high = qc.gaussxw(1024)
# func = lambda x: basis(x) * np.log(np.abs(x - sing_pt))
# exact = 0.0243191054613544
func = lambda x: basis(x) * (x - sing_pt) / ((x - sing_pt) ** 3)
exact = -0.46579
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

n = 10
step = 4.0
hs = 0.25 * ((1.0 / step) ** np.arange(0, n))
print hs

ests = [[]]

total_pts = 16
for m in range(n):
    if m > 0:
        outer = hs[m - 1]
        previous = ests[0][-1]
        # exact = ests[-1][-1]
    else:
        # exact = 0
        outer = 0.5
        previous = with_hole(0.5, gauss[17])
        previous_exact = with_hole(0.5, gauss_high)
        print abs(previous - previous_exact)
    pts = 8 * 24 + 1
    total_pts += pts
    addt = donut(outer, hs[m], gauss[pts])
    addt_ex = donut(outer, hs[m], gauss_high)
    # for hypersingular
    if m > 0:
        print (2 * basis(sing_pt) / hs[m]) - (2 * basis(sing_pt) / hs[m-1])
    print hs[m]
    print addt
    print addt_ex
    assert(abs(addt - addt_ex) < 1e-8)
    current = previous + addt
    ests[0].append(current)
    # ests.append([])
    # for i in range(1, m + 1):
    #     factor = 1.0 / ((step ** (i) - 1.0)
    #     prev_bad = ests[i - 1][-2]
    #     prev_good = ests[i - 1][-1]
    #     new = factor * ((step ** i) * prev_good - prev_bad)
    #     ests[i].append(new)
    # print ests
    best_yet = ests[-1][-1]
    best_no_rich = ests[0][-1]
    error = abs(exact - best_yet)
    no_rich_err = abs(exact - best_no_rich)
    print "After " + str(m) + " iterations: " + str(error) + "     with value: " + str(best_yet)
    print "Without Richardson, best is: " + str(no_rich_err) + "    with value: " + str(best_no_rich)

initial = ests[0]
import ipdb; ipdb.set_trace()
print total_pts
