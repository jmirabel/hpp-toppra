import toppra.cpp
import numpy as np
import matplotlib.pyplot as plt

# P(t) = sum_i=0..5( c_i t^i )
# P(0) = 0 = c_0
# P'(0) = 0 = c_1
# P''(0) = 0 = c_2

# P(1) = 1 = c_3 + c_4 + c_5
# P'(1) = 0 = 3 c_3 + 4 c_4 + 5 c_5
# P''(1) = 0 = 6 c_3 + 12 c_4 + 20 c_5
# c_3 = 10
# c_4 = -15
# c_5 = 6

K = np.zeros((6, 6))
t_order_conditions = [ (0,0) for _ in range(6) ]
boundary_conditions = np.zeros(6)


def set_boundary_condition(row, order, s, value):
    """Set the row-th boundary condition:

    d^order P / ds^order (s) = value
    """
    assert order in [0, 1, 2]
    if order == 0:
        C = [
            1,
        ] * 6
    if order == 1:
        C = [0, 1, 2, 3, 4, 5]
    if order == 2:
        C = [0, 0, 2, 6, 12, 20]

    d = K.shape[1] - 1
    for j in range(order):
        K[row, d - j] = 0
    K[row, d - order] = 1
    for j in range(order + 1, K.shape[1]):
        K[row, d - j] = s * K[row, d - j + 1]

    K[row] *= list(reversed(C))
    t_order_conditions[row] = (s, order)
    boundary_conditions[row] = value


set_boundary_condition(0, 0, 0.0, 0.0)  # P(0) = 0
set_boundary_condition(1, 0, 1.0, 1.0)  # P(1) = 1
set_boundary_condition(2, 1, 0.0, 0.0)  # P'(0) = 0
set_boundary_condition(3, 1, 1.0, 0.0)  # P'(1) = 0
set_boundary_condition(4, 2, 0., 0.0)  # P''(0) = 0
# set_boundary_condition(4, 1, 0.5, 0.0)  # P'(0.5) = 0
set_boundary_condition(5, 2, 1.0, 0.0)  # P''(1) = 0

spline_coeffs = np.linalg.solve(K, boundary_conditions)

path = toppra.cpp.PiecewisePolyPath([spline_coeffs.tolist()], [0, 1])
for (t, order), value in zip(t_order_conditions, boundary_conditions):
    np.testing.assert_allclose(path.eval_single(t, order), value, atol=1e-8)

vel = toppra.cpp.LinearJointVelocity([-1], [1])
vel.discretizationType = toppra.cpp.Interpolation
acc = toppra.cpp.LinearJointAcceleration([-2], [2])
acc.discretizationType = toppra.cpp.Interpolation

algo = toppra.cpp.TOPPRA([vel, acc], path)
algo.setN(200)
algo.computePathParametrization()

# Reconstruct s, sd and sdd
s = algo.parametrizationData.gridpoints
ds = s[1:] - s[:-1]
sd2 = algo.parametrizationData.parametrization
sd = np.sqrt(sd2)
sdd = (sd2[1:] - sd2[:-1]) / 2 / ds

# Calculate t
dt = ds / (sd[1:] + sd[:-1]) * 2
t = np.insert(np.cumsum(dt), 0, 0.0)

# Assuming constant acceleration, calculate spline coefficient
# s_i(t) = a[i] + b[i] (t - t0[i]) + c[i] (t - t0[i])**2
class ConstAccel:
    def __init__(self, t, s, sd, path):
        sd2 = sd**2
        ds = s[1:] - s[:-1]
        self._t0 = t[:-1]
        self._c = (sd2[1:] - sd2[:-1]) / ds / 4
        self._b = sd[:-1]
        self._a = s[:-1]
        self._q = path

    def s(self, t, order=0):
        for i, t0 in enumerate(self._t0):
            if t < t0:
                if i > 0:
                    i = i - 1
                break
        dt = t - self._t0[i]
        if order == 0:
            return self._a[i] + self._b[i] * dt + self._c[i] * dt**2
        elif order == 1:
            return self._b[i] + 2 * self._c[i] * dt
        elif order == 2:
            return 2 * self._c[i]
        else:
            return 0.0

    def q(self, t, order=0):
        assert order in [0, 1, 2]
        s = self.s(t, order=0)
        q = path.eval_single(s, order=0)
        if order == 0:
            return q
        sd = self.s(t, order=1)
        qs = path.eval_single(s, order=1)
        if order == 1:
            return qs * sd
        sdd = self.s(t, order=2)
        qss = path.eval_single(s, order=2)
        if order == 2:
            return qs * sdd + qss * sd**2

ss = np.linspace(0, 1, num=1000)
q = [path.eval_single(u, order=0) for u in ss]
v = [path.eval_single(u, order=1) for u in ss]
a = [path.eval_single(u, order=2) for u in ss]
plt.figure()
plt.title("Initial trajectory")
plt.plot(ss, q, label="q")
plt.plot(ss, v, label="v")
plt.plot(ss, a, label="a")
plt.legend()

S = ConstAccel(t, s, sd, path)

ts = np.linspace(0, t[-1], num=10000)
q = [S.q(u, order=0) for u in ts]
v = [S.q(u, order=1) for u in ts]
a = [S.q(u, order=2) for u in ts]
plt.figure()
plt.title("Timed trajectory")
plt.plot(ts, q, label="q")
plt.plot(ts, v, label="v")
plt.plot(ts, a, label="a")
plt.legend()

plt.figure()
plt.title("Time parametrization")
#plt.plot(ts, [ S.s(u, order=0) for u in ts ], label="s")
plt.plot(ts, [ S.s(u, order=1) for u in ts ], label="sd")
#plt.plot(ts, [ S.s(u, order=2) for u in ts ], label="sdd")
plt.legend()
plt.show()
