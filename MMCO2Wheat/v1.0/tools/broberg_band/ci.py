# 95 % CI of the 420->550 ppm grain-yield response from the rebuilt Broberg et al. (2019) data
import numpy as np
P = np.load("pts.npy"); x, y = P[:, 0], P[:, 1]
def resp(b): f = lambda c: np.polyval(b, c); return 100 * ((1 + f(550) / 100) / (1 + f(420) / 100) - 1)
def drop(x, y, k):
    keep = np.ones(len(x), bool)
    for _ in range(k):
        b = np.polyfit(x[keep], y[keep], 2); r = np.abs(y - np.polyval(b, x)); r[~keep] = -1; keep[np.argmax(r)] = False
    return x[keep], y[keep]
rng = np.random.default_rng(1)
for k in (13, 15):
    xs, ys = drop(x, y, k)
    for deg in (1, 2):
        b, V = np.polyfit(xs, ys, deg, cov=True)
        draws = rng.multivariate_normal(b, V, 20000)
        par = np.percentile([resp(d) for d in draws], [2.5, 97.5])
        boot = []
        for _ in range(5000):
            i = rng.integers(0, len(xs), len(xs)); boot.append(resp(np.polyfit(xs[i], ys[i], deg)))
        print("outliers removed %d, %s: est %.1f %%, 95%% CI parametric %.1f-%.1f, bootstrap %.1f-%.1f"
              % (k, ["", "linear", "quadratic"][deg], resp(b), *par, *np.percentile(boot, [2.5, 97.5])))
# Broberg's published coefficients for reference
print("Broberg Table 1: linear %.1f, quadratic %.1f" % (resp([0.062, -17.5]), resp([-0.00049, 0.576, -142.9])))
