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



def eval(d, q):
    return integrate(K([outer_pt[0], outer_pt[1] + d]), q)

quad_high = qc.gaussxw(5000)
outer_pt = [0.2, 0.00]
# exact = 0.311900553157219
exact = -0.41656692028

for n_q in 2 ** np.arange(1, 13):
    print "Quad pts: " + str(n_q)
    start_dist = 7 * (2.0 / n_q)
    quad_low = qc.gaussxw(n_q)

    min_error = 1000
    min_val = 0
    for n in range(2, 20):
        ds = np.linspace(start_dist, 1.0, n)
        init = [eval(d, quad_low) for d in ds]
        # init_exact = np.array([eval(d, quad_high) for d in ds])
        poly = spi.lagrange(ds, init)
        taylor_exp = []
        result = poly(0.0)
        error = abs(result - exact)
        if min_error > error:
            min_error = error
            min_val = (n, result)
    print "Min at order expansion: " + str(min_val[0])
    print "Min error is : " + str(min_error)
    print "Value is : " + str(min_val[1])
    print "\n\n"

