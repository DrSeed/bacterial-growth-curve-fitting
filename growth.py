# Bacterial growth curves: model-free max slope and a Gompertz fit (Zwietering form).
import numpy as np


def gompertz(t, A, mu, lam):
    # y = ln(N/N0): asymptote A, max specific growth rate mu, lag time lam.
    return A * np.exp(-np.exp(mu * np.e / A * (lam - t) + 1))


def log_od(od, blank=0.0):
    return np.log(np.clip(np.asarray(od) - blank, 1e-4, None))


def max_slope(t, od, blank=0.0, window=8):
    """Sliding-window regression on ln(OD); returns (mu_max, lag, doubling_time)."""
    y = log_od(od, blank)
    best = (-np.inf, 0.0, 0.0)
    for i in range(len(t) - window + 1):
        tt, yy = t[i:i + window], y[i:i + window]
        slope, intercept = np.polyfit(tt, yy, 1)
        if slope > best[0]:
            best = (slope, tt.mean(), slope * tt.mean() + intercept)
    mu, t_star, y_star = best
    lag = t_star - (y_star - y[:3].mean()) / mu
    return mu, lag, np.log(2) / mu


def fit_gompertz(t, od, blank=0.0, iters=300):
    """Levenberg-Marquardt fit of the Gompertz model to ln(OD/OD0)."""
    y = log_od(od, blank)
    y = y - y[:3].mean()
    mu0, lag0, _ = max_slope(t, od, blank)
    p = np.array([max(y.max(), 0.5), mu0, max(lag0, 0.0)])
    lam_lm = 1e-2

    def sse(q):
        return float(np.sum((y - gompertz(t, *q)) ** 2))

    cur = sse(p)
    for _ in range(iters):
        J = np.empty((len(t), 3))
        for k in range(3):
            h = 1e-6 * max(abs(p[k]), 1e-3)
            dp = p.copy()
            dp[k] += h
            J[:, k] = (gompertz(t, *dp) - gompertz(t, *p)) / h
        r = y - gompertz(t, *p)
        A = J.T @ J
        step = np.linalg.solve(A + lam_lm * np.diag(np.diag(A) + 1e-12), J.T @ r)
        new = p + step
        if new[0] > 0 and new[1] > 0 and sse(new) < cur:
            p, cur = new, sse(new)
            lam_lm *= 0.3
        else:
            lam_lm *= 10
        if lam_lm > 1e12:
            break
    A_fit, mu, lag = p
    return {"A": A_fit, "mu": mu, "lag": lag, "doubling_time": np.log(2) / mu}
