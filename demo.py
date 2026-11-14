# Self-contained demo: three strains on a simulated plate reader, analysed with and without blank subtraction.
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from growth import gompertz, max_slope, fit_gompertz, log_od

BLANK = 0.045
STRAINS = {"wild type": (0.90, 1.5), "slow mutant": (0.50, 1.0), "long-lag mutant": (0.90, 4.0)}


def plate_reader(mu, lam, rng, od0=0.01, hours=16, step=0.25):
    t = np.arange(0, hours + 1e-9, step)
    A = np.log(1.2 / od0)
    od = od0 * np.exp(gompertz(t, A, mu, lam))
    return t, od * rng.lognormal(0, 0.02, len(t)) + BLANK + rng.normal(0, 0.002, len(t))


def main():
    os.makedirs("figures", exist_ok=True)
    os.makedirs("results", exist_ok=True)
    rng = np.random.default_rng(21)
    fig, axes = plt.subplots(1, 3, figsize=(17, 4.8))
    rows = []
    for (name, (mu, lam)), c in zip(STRAINS.items(), ["#1f77b4", "#ff7f0e", "#2ca02c"]):
        t, od = plate_reader(mu, lam, rng)
        g = fit_gompertz(t, od, blank=BLANK)
        ms = max_slope(t, od, blank=BLANK)
        raw = max_slope(t, od, blank=0.0)
        rows.append({"strain": name, "true_mu": mu, "true_lag_h": lam,
                     "gompertz_mu": round(g["mu"], 3), "gompertz_lag_h": round(g["lag"], 2),
                     "max_slope_mu": round(ms[0], 3), "max_slope_lag_h": round(ms[1], 2),
                     "doubling_time_min": round(60 * g["doubling_time"], 1),
                     "no_blank_mu": round(raw[0], 3), "no_blank_lag_h": round(raw[1], 2)})
        axes[0].plot(t, od, ".", color=c, ms=3, label=name)
        y = log_od(od, BLANK)
        axes[1].plot(t, y - y[:3].mean(), ".", color=c, ms=3)
        axes[1].plot(t, gompertz(t, g["A"], g["mu"], g["lag"]), color=c,
                     label=f"{name}: mu={g['mu']:.2f}/h, lag={g['lag']:.1f} h")
    table = pd.DataFrame(rows)
    table.to_csv("results/summary.csv", index=False)
    axes[0].axhline(BLANK, color="k", ls=":", label="medium blank")
    axes[0].set_xlabel("time (h)")
    axes[0].set_ylabel("OD600 (raw)")
    axes[0].set_title("Raw plate-reader curves")
    axes[0].legend(fontsize=8)
    axes[1].set_xlabel("time (h)")
    axes[1].set_ylabel("ln(OD/OD0), blank subtracted")
    axes[1].set_title("Gompertz fits")
    axes[1].legend(fontsize=8)

    x = np.arange(len(table))
    axes[2].bar(x - 0.27, table["true_mu"], 0.27, label="truth", color="#444444")
    axes[2].bar(x, table["gompertz_mu"], 0.27, label="blank subtracted", color="#2a7f62")
    axes[2].bar(x + 0.27, table["no_blank_mu"], 0.27, label="no blank subtraction", color="#c0504d")
    axes[2].set_xticks(x, table["strain"])
    axes[2].set_ylabel("max growth rate mu (per h)")
    axes[2].set_title("Skipping the blank flattens growth")
    axes[2].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig("figures/demo.png", dpi=120)
    print(table.to_string(index=False))


if __name__ == "__main__":
    main()
