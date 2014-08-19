import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import quadracheer as qc
from scipy.special import hankel1, jn
from scipy.special import legendre
import scipy.interpolate as spi
from functools import partial
from math import factorial
matplotlib.rcParams['font.family'] = 'serif'
matplotlib.rcParams['font.serif'] = ['Computer Modern Roman']
matplotlib.rcParams['text.usetex'] = True
matplotlib.rcParams['font.size'] = 14
matplotlib.rcParams['xtick.direction'] = 'out'
matplotlib.rcParams['ytick.direction'] = 'out'
matplotlib.rcParams['lines.linewidth'] = 3

def integrate(f, quad_rule):
    sum = 0
    for (p, w) in zip(quad_rule[0], quad_rule[1]):
        sum += w * f([p, 0.0])
    return sum

def dist(x1,x2): return np.sqrt((x2[0] - x1[0]) ** 2 + (x2[1] - x1[1]) ** 2)
# The laplace single layer potential
single_layer = lambda x1: lambda x2: (-1.0 / (2 * np.pi)) * np.log(dist(x1,x2))

# The laplace double layer potential
double_layer = lambda x1: lambda x2: (x2[0] - x1[0]) / (dist(x1, x2) ** 2)

# Hypersingular thingamabob
hypersing = lambda x1: lambda x2: (x2[0] - x1[0]) / (dist(x1, x2) ** 3)

quad_high = qc.gaussxw(2000)
def test(n_q, n, outer_pt, K, exact, step = 2.0, dist_mult = 5):
    quad_low = qc.gaussxw(n_q)
    start = dist_mult * (2.0 / n_q)
    ds = (start * (step ** np.arange(0, n))).tolist()
    ds.reverse()

    def eval(d, q):
        return integrate(K([outer_pt[0], outer_pt[1] + d]), q)

    init = np.array([eval(d, quad_low) for d in ds])
    # print(str(init))
    rich = [init]
    for m in range(1, n):
        prev_rich = rich[-1]
        best_est = prev_rich[-1]
        error = abs(best_est - exact)
        # print("Tier " + str(m - 1) + " equals: " + str(best_est) +
        #         "   with an error of: " + str(error))
        factor = (1.0 / (step ** m - 1.0))
        next_rich = factor * ((step ** m) * prev_rich[1:] - prev_rich[:-1])
        rich.append(next_rich)
    print rich

    final_est = rich[-1][-1]
    error = abs(final_est - exact)
    error_estimate = 2 * abs(final_est - rich[-2][-1])
    return error

# n = 9
# n_q = 128
# outer_pt = [0.2, 0.0]
# exact = 0.311900553157219
# kernel = single_layer
# outer_pt = [-0.5, 0.0]
# exact = 0.276671121911897
# kernel = single_layer
# outer_pt = [0.5, 0.0]
# exact = -1.09861228866811
# kernel = double_layer
outer_pt = [0.02, 0.0]
exact = -0.294837925881334
basis = legendre(40)
kernel = lambda x1: lambda x2: basis(x2[0]) * double_layer(x1)(x2)
# for n_q in 2 ** np.arange(1, 10):
#     print test(n_q, n, single_layer, exact)

ns = [12]
n_q = 240
steps = [1.2]
dms = [7.0]
best = 1000
for dist_mult in dms:
    # for n in range(6, 16, 2):
    for n in ns:
        for step in steps:
            e = test(n_q, n, outer_pt, kernel,
                     exact, step = step, dist_mult = dist_mult)
            if e < best:
                best = e
                best_info = (step, n, dist_mult)
print best
print best_info
