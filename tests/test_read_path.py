import toppra.cpp
import toppra.interpolator
import hpp_toppra_cpp
import matplotlib.pyplot as plt
import numpy as np
import typing as T

from eureka.config.device_config.robots import URRobotConfig

robot_config = URRobotConfig(model="3e")
vel_limits = robot_config.joint_velocity_limits[1, :]
acc_limits = robot_config.joint_acceleration_limits[1, :]

device = hpp_toppra_cpp.pinocchio.Device.create("world")
hpp_toppra_cpp.pinocchio.urdf.loadModel(
    device,
    0,
    "",
    "anchor",
    "package://ur_description_eureka/urdf/ur3e.urdf",
    "package://ur_description_eureka/srdf/ur3e.srdf",
)

path = hpp_toppra_cpp.core.readBinPath(device, "../last-hpp-toppra-path.bin")
print(path)

times = [0]
subpaths = []
for i in range(path.numberPaths()):
    subpath = path.pathAtRank(i)
    subpaths.append(subpath)
    times.append(times[-1] + subpath.timeRange[1])
n = len(path.eval(0)[0])


def draw_vertical_lines(times):
    for t in times:
        plt.axvline(t)


def plot_path(path):
    times = []
    if isinstance(path, hpp_toppra_cpp.core.PathVector):
        times.append(0)
        for i in range(path.numberPaths()):
            subpath = path.pathAtRank(i)
            times.append(times[-1] + subpath.timeRange[1])

    ts = np.linspace(path.timeRange[0], path.timeRange[1], num=10000)
    plt.subplot(3, 1, 1)
    lines = plt.plot(ts, [path.eval(t)[0] for t in ts])
    for i, l in enumerate(lines):
        l.set_label(f"q{i}")
    draw_vertical_lines(times)
    plt.legend()
    plt.subplot(3, 1, 2)
    lines = plt.plot(ts, [path.derivative(t, 1) for t in ts])
    for i, l in enumerate(lines):
        l.set_label(f"v{i}")
    draw_vertical_lines(times)
    plt.legend()
    plt.subplot(3, 1, 3)
    lines = plt.plot(ts, [path.derivative(t, 2) for t in ts], "+")
    for i, l in enumerate(lines):
        l.set_label(f"a{i}")
    draw_vertical_lines(times)
    plt.legend()


# Assuming constant acceleration, calculate spline coefficient
# s_i(t) = a[i] + b[i] (t - t0[i]) + c[i] (t - t0[i])**2
class ConstAccel:
    def __init__(self, t, s, sd, path: toppra.cpp.GeometricPath):
        sd2 = sd**2
        ds = s[1:] - s[:-1]
        self._t0 = t[:-1]
        self._T = t[-1]
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
        q = self._q.eval_single(s, order=0)
        if order == 0:
            return q
        sd = self.s(t, order=1)
        qs = self._q.eval_single(s, order=1)
        if order == 1:
            return qs * sd
        sdd = self.s(t, order=2)
        qss = self._q.eval_single(s, order=2)
        if order == 2:
            return qs * sdd + qss * sd**2

    def eval(self, t):
        return self.q(t, order=0), True

    def derivative(self, t, order):
        return self.q(t, order=order)

    @property
    def timeRange(self):
        return 0, self._T


def run_toppra(
    path: hpp_toppra_cpp.core.Path,
    vel_limits: T.Optional[np.ndarray] = None,
    acc_limits: T.Optional[np.ndarray] = None,
    N: int = 200,
):
    path_wrapper = hpp_toppra_cpp.core.TOPPRAPathWrapper(path)
    constraints = []
    if vel_limits is not None:
        vel = toppra.cpp.LinearJointVelocity(-vel_limits, vel_limits)
        vel.discretizationType = toppra.cpp.Interpolation
        constraints.append(vel)
    if acc_limits is not None:
        acc = toppra.cpp.LinearJointAcceleration(-acc_limits, acc_limits)
        acc.discretizationType = toppra.cpp.Interpolation
        constraints.append(acc)

    algo = toppra.cpp.TOPPRA(constraints, path_wrapper)
    gridpointsMethod = 2
    if gridpointsMethod == 0:
        algo.setGridpoints(toppra.interpolator.propose_gridpoints(path_wrapper,
            min_nb_points=N))
    elif gridpointsMethod == 1:
        algo.setGridpoints(path_wrapper.proposeGridpoints(minNbPoints=N))
    elif gridpointsMethod == 2:
        initialGridpoints = [0.]
        for i in range(path.numberPaths()):
            initialGridpoints.append(initialGridpoints[-1] + path.pathAtRank(i).length())
        algo.setGridpoints(path_wrapper.proposeGridpoints(minNbPoints=N, initialGridpoints=initialGridpoints))
    else:
        algo.setN(N)
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

    S = ConstAccel(t, s, sd, path_wrapper)

    return S, (t, s, sd, sdd)


# plot_path(path)
# plt.show()

S, rets = run_toppra(path, vel_limits, acc_limits, N=1000)
plot_path(S)
plt.show()
