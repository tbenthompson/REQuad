import matplotlib.pyplot as plt
import numpy as np
import quadracheer as qc
from scipy.special import legendre

def dist(x1,x2): return np.sqrt((x2[0] - x1[0]) ** 2 + (x2[1] - x1[1]) ** 2)
single_layer = lambda x1: lambda x2: (-1.0 / (2 * np.pi)) * np.log(dist(x1,x2))
double_layer = lambda x1: lambda x2: (x2[0] - x1[0]) / (dist(x1, x2) ** 2)
hypersing = lambda x1: lambda x2: (x2[0] - x1[0]) * (x1[1] - x2[1]) / (dist(x1, x2) ** 4)

basis = legendre(16)
kernel = lambda x1: lambda x2: basis(x2[0]) * single_layer(x1)(x2)
# kernel = lambda x1: lambda x2: basis(x2[0]) * double_layer(x1)(x2)
# kernel = lambda x1: lambda x2: basis(x2[0]) * hypersing(x1)(x2)

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


def richardson(values, step, log = False):
    rich = [values]
    for m in range(1, len(values)):
        # if m == 1 and log is True:
        #     factor =
        prev_rich = rich[m - 1]
        mult = step ** (2 * m - 1)
        factor = (1.0 / (mult - 1.0))
        next_rich = factor * (mult * prev_rich[1:] - prev_rich[:-1])
        rich.append(next_rich)
    return rich


def simple():

    high_quad_order = 500
    highquad = qc.gaussxw(high_quad_order)
    med_quad_order = 100
    medquad = qc.gaussxw(med_quad_order)

    low_quad_order = 16
    lowquad = qc.gaussxw(low_quad_order)
    safe_dist = 6.0 / low_quad_order
    print safe_dist

    n = 8
    step = 2.0
    hs = safe_dist * ((1.0 / step) ** np.arange(0, n))
    print hs

    sing_pt = 0.2
    d = 0.01
    perfect = -0.00580747813511577
    # perfect = 0.113028326740308
    # perfect = 0.505981109454

    integrals = []
    for i, h in enumerate(hs):
        f = kernel([sing_pt, h])
        est = recursive(sing_pt, hs[:max(2, i)], f, lowquad)
        exact = recursive(sing_pt, hs[:max(2, i)], f, medquad)
        error = abs(est - exact)
        print error
        integrals.append(est)
    integrals = np.array(integrals)
    # print integrals[-1]
    print integrals[-1] - perfect

    rich = richardson(integrals, step)
    # print rich[-1][-1]
    print rich[-3][-1] - perfect
    plt.loglog(hs, np.abs(integrals - perfect), label = 0)
    for i in range(1, len(rich)):
        plt.loglog(hs[i:], np.abs(rich[i] - perfect), '.-', linewidth = 2, label = i)
    plt.legend()
    plt.show()

def rich_fun():
    src_pt = [0.0, -0.5]
    K = kernel(src_pt)
    I = K([0.0, 0.0])
    hs = np.linspace(0.000, 1.0, 101.0)[1:]
    I_est = np.array([K([0.0, h]) for h in hs])
    hs_ratio = hs[1:] / hs[:-1]
    print hs_ratio
    better_est = (hs_ratio * I_est[:-1] - I_est[1:]) / (hs_ratio - 1)
    even_better_est = (hs_ratio[:-1] ** 2 * better_est[:-1] - better_est[1:]) / (hs_ratio[:-1] ** 2 - 1)
    import ipdb; ipdb.set_trace()
    plt.plot(hs, I_est - I, 'y-')
    plt.plot(hs[:-1], better_est - I, 'g-')
    plt.plot(hs[:-2], even_better_est - I, 'b-')
    plt.show()

if __name__ == '__main__':
    simple()
    # rich_fun()
