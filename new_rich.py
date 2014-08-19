import matplotlib.pyplot as plt
import numpy as np
from scipy.special import legendre
from gaussian_quad import gaussxw
from scipy.interpolate import BarycentricInterpolator as bi
import sympy as sp



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

def plot_best_rich(hs, rich, perfect, which):
    best = [abs(r[0] - perfect) for r in rich]
    fig, ax = plt.subplots(1)
    ax.loglog(hs, best, 'o-', linewidth = 2, label = "RICH")
    ax.set_xlabel('$h$')
    ax.set_ylabel('$E$')
    fig.savefig('error_' + which + '.pdf')

def run(problem, name):
    K, sing_pt, basis, perfect, include_pt = problem

    high_quad_order = 500
    highquad = gaussxw(high_quad_order)
    med_quad_order = 100
    medquad = gaussxw(med_quad_order)
    low_quad_order = 16
    lowquad = gaussxw(low_quad_order)

    safe_dist = 6.0 / low_quad_order

    n = 4
    step = 2.0
    hs = safe_dist * ((1.0 / step) ** np.arange(0, n))
    print("Distances are: " + str(hs))

    print("Correct answer is: " + str(perfect))
    kernel = lambda x1: lambda x2: basis(x2[0]) * K(x1[0], x2[0], x1[1], x2[1])

    integrals = []
    for i, h in enumerate(hs):
        f = kernel([sing_pt, h])
        est = recursive(sing_pt, hs[:], f,
                        lowquad, include_pt = include_pt)
        exact = recursive(sing_pt, hs[:], f,
                          medquad, include_pt = include_pt)
        error = abs(est - exact)
        print("Quadrature error " + str(i) + ": " + str(error) + "   and " +
              "Distance error " + str(i) + ": " + str(exact - perfect))
        integrals.append(est)
    integrals = np.array(integrals)
    no_rich_error = integrals[-1] - perfect
    print("Raw integral error: " + str(no_rich_error))

    rich = richardson(hs, integrals, step)
    print("Best richardson error: " + str(abs(rich[-1][-1] - perfect))) + "\n\n"

    interp_est = interp(integrals, hs)
    # print("Interpolation error: " + str(np.array(interp_est) - perfect))

    # plot_all_rich(hs, rich, perfect)
    plot_best_rich(hs, rich, perfect, name)

    # This prints a random number
    print 88

def main():
    """
    Start with the poisson equation single layer potential
    Take normal derivative with respect to observation variable to get double
    layer potential.
    Take normal derivative with respect to source/integration variable to get
    hypersingular potential
    Normal derivatives are in the y direction because the element is along
    the x-axis
    """
    x1, x2, y1, y2 = sp.symbols('x1, x2, y1, y2')
    dist = sp.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2)
    slp = -1 / (2 * sp.pi) * sp.log(dist)
    dlp = sp.diff(slp, y1)
    hlp = sp.diff(dlp, y2)
    args = (x1, x2, y1, y2)
    single_layer = sp.utilities.lambdify(args, slp)
    double_layer = sp.utilities.lambdify(args, dlp)
    hypersing = sp.utilities.lambdify(args, hlp)

    test_problems = dict()
    # Problem format: (kernel, singular_pt, basis, exact, include_pt)
    # include_pt indicates whether to include the nearest point on the element
    # in the integration. For some highly singular integrals, ignoring this
    # point does not hurt convergence and is much more numerically stable.
    test_problems['single1'] = (single_layer, 0.2, legendre(1), 0.0628062411975970, True)
    test_problems['single16'] = (single_layer, 0.2, legendre(16), -0.00580747813511577, True)
    test_problems['double1'] = (double_layer, 0.2, legendre(1), 0.0, True)
    test_problems['double16'] = (double_layer, 0.2, legendre(16), 0.0, True)
    test_problems['hyper1'] = (hypersing, 0.2, legendre(1), -0.130846, True)
    test_problems['hyper16'] = (hypersing, 0.2, legendre(16), 0.0, True)
    # run(test_problems['single1'], 'single1')
    # run(test_problems['single16'], 'single16')
    # run(test_problems['double1'], 'double1')
    # run(test_problems['double16'], 'double16')
    run(test_problems['hyper1'], 'hyper1')
    run(test_problems['hyper16'], 'hyper16')

if __name__ == '__main__':
    main()
