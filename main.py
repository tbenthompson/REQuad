import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import quadracheer as qc
from scipy.special import hankel1, jn
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

def qbx_expand(kernel, quad_rule, exp_pt, order):
    c = []
    h = 1e-3
    left_pt = [exp_pt[0], exp_pt[1] - h]
    lefter_pt = [exp_pt[0], exp_pt[1] - h - h]
    right_pt = [exp_pt[0], exp_pt[1] + h]
    righter_pt = [exp_pt[0], exp_pt[1] + h + h]
    K_c = integrate(kernel(exp_pt), quad_rule)
    K_l = integrate(kernel(right_pt), quad_rule)
    K_r = integrate(kernel(left_pt), quad_rule)
    K_ll = integrate(kernel(righter_pt), quad_rule)
    K_rr = integrate(kernel(lefter_pt), quad_rule)
    const = K_c
    deriv = (-K_rr + 8 * K_r - 8 * K_l + K_ll) / (12 * h)
    second_deriv = (-K_rr + 16 * K_r - 30 * K_c + 16 * K_l - K_ll) / (12 * h ** 2)
    third_deriv = (K_rr - 2 * K_r + 2 * K_l - K_ll) / (2 * h ** 3)
    fourth_deriv = (K_rr - 4 * K_r + 6 * K_c - 4 * K_l + K_ll) / (h ** 4)
    c = [const, deriv, second_deriv, third_deriv, fourth_deriv]
    print c
    return (c, exp_pt)

def qbx_eval(pt, qbx_info):
    return sum([(cm / factorial(m)) * (qbx_info[1][1] - pt[1]) ** m
                for (m, cm) in enumerate(qbx_info[0])])

qbx_quad_order = 16
qbx_quad_factor = 4
quad_orders = [qbx_quad_order, qbx_quad_order * qbx_quad_factor * 4, 5000]
quad_rules = [(qc.gaussxw(n), str(n) + " order gauss") for n in quad_orders]

def dist(x1,x2): return np.sqrt((x2[0] - x1[0]) ** 2 + (x2[1] - x1[1]) ** 2)
# The laplace single layer potential
# K = lambda x1: lambda x2: (-1.0 / (2 * np.pi)) * np.log(dist(x1,x2))
# The laplace double layer potential
K = lambda x1: lambda x2: (x2[0] - x1[0]) / (dist(x1, x2) ** 3)


# Eval on a nice line!
interval_length = 2.0
qbx_distance = 5 * (interval_length / qbx_quad_order)
print qbx_distance
exp_pt = [0.2, qbx_distance]
qbx_info = qbx_expand(K, qc.gaussxw(qbx_quad_order * qbx_quad_factor), exp_pt, 0)
view = [0.005, 0.05]
pts = np.linspace(view[0], view[1], 50)
results = [[integrate(K([exp_pt[0], p]), q[0]) for p in pts] for q in quad_rules]
qbx_res = [qbx_eval([exp_pt[0], p], qbx_info) for p in pts]
quad_rules.insert(-1, (0, "QBX"))
results.insert(-1, qbx_res)
colors = ['r-', 'b-', 'k-']
for i, (r, q) in enumerate(zip(results, quad_rules)[:-1]):
    error = np.abs(np.array(r) - np.array(results[-1]))
    plt.plot(pts, (np.log(error / np.abs(np.array(r))) / np.log(10)), colors[i], label=q[1], linewidth = 3)
plt.legend(loc = 'lower left')
plt.xlabel('$x$')
plt.ylabel('$\log_{10}(e)$')
plt.savefig('qbx_vs_gauss.pdf')
plt.show()

# np.testing.assert_almost_equal(qbx_eval(exp_pt, qbx_info), integrate(K(exp_pt), quad_rules[1][0]), 3)

# Evaluate the integral for a big grid of observation points.
# nx = 100
# ny = 100
# view = [-6.0, 6.0, 6.0, -6.0]
# x_outer = np.linspace(view[0], view[1], nx)
# y_outer = np.linspace(view[2], view[3], ny)
# X, Y = np.meshgrid(x_outer, y_outer)
# F_est = np.zeros_like(X)
# F_ex = np.zeros_like(X)
# error_plot = True
# for i in range(nx):
#     for j in range(ny):
#         pt = [X[i, j], Y[i, j]]
#         if np.sqrt((pt[0] - exp_pt[0]) ** 2 + (pt[1] - exp_pt[1]) ** 2) < exp_pt[1]:
#             F_est[i, j] = qbx_eval(pt, qbx_info)
#         else:
#             F_est[i, j] = integrate(K(pt), q_low)
#         if error_plot:
#             F_ex[i, j] = integrate(K(pt), q_high)
#
# plt.figure()
# if error_plot:
#     plt.imshow(np.log(np.abs(F_est - F_ex)) / np.log(10))
# else:
#     plt.imshow(F_est)
# plt.colorbar()
# plt.show()
