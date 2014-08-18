import matplotlib.pyplot as plt
import numpy as np
import quadracheer as qc
from scipy.special import legendre
from scipy.interpolate import BarycentricInterpolator as bi

def dist(x1,x2): return np.sqrt((x2[0] - x1[0]) ** 2 + (x2[1] - x1[1]) ** 2)
single_layer = lambda x1: lambda x2: (-1.0 / (2 * np.pi)) * np.log(dist(x1,x2))
double_layer = lambda x1: lambda x2: (x2[0] - x1[0]) / (dist(x1, x2) ** 2)
superhyper = lambda x1: lambda x2: (x2[0] - x1[0]) / (dist(x1, x2) ** 8)
hypersing = lambda x1: lambda x2: (x2[0] - x1[0]) / (dist(x1, x2) ** 4)

test_problems = dict()
test_problems['single1'] = (single_layer, legendre(1), 0.0628062411975970, True)
test_problems['single16'] = (single_layer, legendre(16), -0.00580747813511577, True)
test_problems['double1'] = (double_layer, legendre(1), 1.91890697837663, True)
test_problems['double16'] = (double_layer, legendre(16), 0.113028326740308, True)
test_problems['super0'] = (superhyper, legendre(0), -0.579967, False)
test_problems['hyper16'] = (hypersing, legendre(0), 0.505981109454, False)

def integrate(x_min, x_max, f, q):
    x = q[0]
    w = q[1]
    new_w = w * ((x_max - x_min) / 2.0)
    new_x = x_min + (x_max - x_min) * ((x + 1.0) / 2.0)
    return sum([new_w[i] * f([new_x[i], 0.0]) for i in range(len(new_w))])

def with_hole(sing_pt, h, f, q):
    return integrate(-1.0, sing_pt - h, f, q) + integrate(sing_pt + h, 1.0, f, q)

def donut(sing_pt, ho, hi, f, q):
    return integrate(sing_pt - ho, sing_pt - hi, f, q) + integrate(sing_pt + hi, sing_pt + ho, f, q)

def recursive(sing_pt, sizes, f, q, include_pt = False):
    def inner(remaining):
        term = donut(sing_pt, remaining[0], remaining[1], f, q)
        if len(remaining) == 2:
            if include_pt:
                term += donut(sing_pt, remaining[1], 0.0, f, q)
            return term
        return term + inner(remaining[1:])
    start = with_hole(sing_pt, sizes[0], f, q)
    rem = inner(sizes)
    return start + rem


def richardson(hs, values, step):
    rich = [values]
    hs_ratio = hs[:-1] / hs[1:]
    for m in range(1, len(values)):
        prev_rich = rich[m - 1]
        mult = hs_ratio[(m - 1):] ** m
        factor = (1.0 / (mult - 1.0))
        next_rich = factor * (mult * prev_rich[1:] - prev_rich[:-1])
        rich.append(next_rich)
    return rich

def interp(values, xs):
    res = []
    for i in range(2, len(xs) + 1):
        L = bi(xs[:i])
        L.set_yi(values[:i])
        res.append(L(0.0))
    return res

def plot_all_rich(hs, rich, perfect):
    for i in range(0, len(rich)):
        plt.loglog(hs[i:], np.abs(rich[i] - perfect), '.-', linewidth = 2, label = i)
    plt.legend()
    plt.show()

def plot_best_rich(hs, rich, perfect):
    best = [abs(r[0] - perfect) for r in rich]
    print hs
    print best
    plt.loglog(hs, best, 'o-', linewidth = 2, label = "RICH")
    plt.xlabel('$h$')
    plt.ylabel('$E$')
    plt.show()

def simple():
    high_quad_order = 500
    highquad = qc.gaussxw(high_quad_order)
    med_quad_order = 100
    medquad = qc.gaussxw(med_quad_order)
    low_quad_order = 8
    lowquad = qc.gaussxw(low_quad_order)

    safe_dist = 6.0 / low_quad_order

    n = 10
    step = 2.0
    hs = safe_dist * ((1.0 / step) ** np.arange(0, n))
    print("Distances are: " + str(hs))

    sing_pt = 0.2
    d = 0.01

    which = 'double1'
    which = 'double16'
    # which = 'super0'
    # which = 'single16'
    # which = 'hyper16'
    K, basis, perfect, include_pt = test_problems[which]
    print("Correct answer is: " + str(perfect))
    kernel = lambda x1: lambda x2: basis(x2[0]) * K(x1)(x2)

    integrals = []
    for i, h in enumerate(hs):
        f = kernel([sing_pt, h])
        est = recursive(sing_pt, hs[:], f,
                        lowquad, include_pt = include_pt)
        exact = recursive(sing_pt, hs[:], f,
                          medquad, include_pt = include_pt)
        error = abs(est - exact)
        print("Integration error " + str(i) + ": " + str(error))
        integrals.append(est)
    integrals = np.array(integrals)
    no_rich_error = integrals[-1] - perfect
    print("Raw integral error: " + str(no_rich_error))

    rich = richardson(hs, integrals, step)
    print("Best richardson error: " + str(rich[-3][-1] - perfect))

    interp_est = interp(integrals, hs)
    # print("Interpolation error: " + str(np.array(interp_est) - perfect))

    # plot_all_rich(hs, rich, perfect)
    plot_best_rich(hs, rich, perfect)

def rich_fun():
    src_pt = [0.0, -0.5]
    K = single_layer(src_pt)
    I = K([0.0, 0.0])
    hs = np.linspace(0.000, 1.0, 101.0)[1:]
    I_est = np.array([K([0.0, h]) for h in hs])
    hs_ratio = hs[1:] / hs[:-1]
    print hs_ratio
    new_ests = [I_est]
    n = 9
    for i in range(n):
        factor = hs_ratio[i:] ** (i + 1)
        better = (factor * new_ests[i][:-1] - new_ests[i][1:]) / (factor - 1)
        new_ests.append(better)
    for i in range(n):
        x = hs[:-i]
        if i is 0:
            x = hs
        plt.plot(x, new_ests[i] - I, label = i, linewidth = 3)
    plt.legend()
    plt.show()

if __name__ == '__main__':
    simple()
    # rich_fun()
