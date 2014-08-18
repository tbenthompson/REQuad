import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import quadracheer as qc
from scipy.special import hankel1, jn
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

quad_orders = [16, 256, 5000]
quad_rules = [(qc.gaussxw(n), str(n) + " order gauss") for n in quad_orders]

def dist(x1,x2): return np.sqrt((x2[0] - x1[0]) ** 2 + (x2[1] - x1[1]) ** 2)
# The laplace single layer potential
# K = lambda x1: lambda x2: (-1.0 / (2 * np.pi)) * np.log(dist(x1,x2))
# The laplace double layer potential
K = lambda x1: lambda x2: (x2[0] - x1[0]) / (dist(x1, x2) ** 3)


n = 9
n_q = 128
outer_pt = [0.2, 0.0]
for n_q in 2 ** np.arange(1, 10):
    quad_low = qc.gaussxw(n_q)
    quad_high = qc.gaussxw(5000)
    # exact = integrate(K([outer_pt[0], outer_pt[1]]), quad_high)
    # exact = 0.311900553157219
    exact = -0.41668101144
    start = 5 * (2.0 / n_q)
    step = 2.0
    ds = (start * (step ** np.arange(0, n))).tolist()
    ds.reverse()
    # print ds

    def eval(d, q):
        return integrate(K([outer_pt[0], outer_pt[1] + d]), q)

    init = np.array([eval(d, quad_low) for d in ds])
    init_exact = np.array([eval(d, quad_high) for d in ds])
    # print np.abs(init - init_exact)
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

    final_est = rich[-1][-1]
    error = abs(final_est - exact)
    error_estimate = 2 * abs(final_est - rich[-2][-1])
    print error
# print("Final estimate: " + str(final_est))
# print("Final error estimate:" + str(error_estimate))
# print("Final error: " + str(error))
# print("Exact: " + str(exact))


# integrate(outer_pt

# view = [0.005, 5.0]
# pts = np.linspace(view[0], view[1], 250)
# results = [[integrate(K([outer_pt[0], p]), q[0]) for p in pts] for q in quad_rules]
#
# colors = ['r-', 'b-', 'k-']
# for i, (r, q) in enumerate(zip(results, quad_rules)[:-1]):
#     error = np.abs(np.array(r) - np.array(results[-1]))
#     plt.plot(pts, (np.log(error / np.abs(np.array(r))) / np.log(10)), colors[i], label=q[1], linewidth = 3)
# plt.legend(loc = 'lower left')
# plt.xlabel('$x$')
# plt.ylabel('$\log_{10}(e)$')
# plt.show()
