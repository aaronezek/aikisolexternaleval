"""Fig 1 (ROC + PR with bootstrap bands) and Fig 2 (reliability) from results/eval_scores_all_models.csv."""
import sys, os, json
sys.path.insert(0, os.path.dirname(__file__))
import numpy as np, pandas as pd, lib
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, precision_recall_curve
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8, "axes.labelsize": 8, "xtick.labelsize": 6.5, "ytick.labelsize": 6.5,
    "legend.fontsize": 7, "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6, "xtick.major.size": 2.5, "ytick.major.size": 2.5,
    "axes.spines.top": False, "axes.spines.right": False, "savefig.dpi": 300, "figure.dpi": 150})
D = pd.read_csv(lib.RES/"eval_scores_all_models.csv"); y = D.label.values
FIG = lib.ROOT/"figures"; FIG.mkdir(exist_ok=True)
MODELS = [("aikisol_marginal", "Aiki-Sol (marginal output)", "#0072B2", 1.8), ("logreg5", "Logistic regression, 5 features", "#E69F00", 1.2), ("swi", "SWI (published)", "#8c8c8c", 1.2)]
B = 1000; idx = lib.boot_indices(len(y), B, lib.SEED); gx = np.linspace(0, 1, 201)
def band(curve_fn, s):
    ys = []
    for ix in idx:
        yy = y[ix]
        if yy.min() == yy.max(): continue
        ys.append(curve_fn(yy, s[ix]))
    ys = np.array(ys); return np.quantile(ys, 0.025, axis=0), np.quantile(ys, 0.975, axis=0)
def roc_on_grid(yy, ss):
    f, t, _ = roc_curve(yy, ss); return np.interp(gx, f, t)
def pr_on_grid(yy, ss):
    p, r, _ = precision_recall_curve(yy, ss); o = np.argsort(r); return np.interp(gx, r[o], p[o])   # precision at recall grid
# ---------------- Figure 1
fig, axs = plt.subplots(1, 2, figsize=(6.8, 3.0))
ax = axs[0]
for k, lab, c, lw in MODELS:
    s = D[k].values; f, t, _ = roc_curve(y, s); a = lib.auroc(y, s); lo, hi = band(roc_on_grid, s)
    ax.fill_between(gx, lo, hi, color=c, alpha=0.18, lw=0); ax.plot(f, t, color=c, lw=lw, label=f"{lab}  ({a:.2f})")
ax.plot([0, 1], [0, 1], ls=":", color="#b0b0b0", lw=0.8); ax.set_xlabel("False positive rate"); ax.set_ylabel("True positive rate")
ax.set_xlim(-0.02, 1.02); ax.set_ylim(-0.02, 1.02); ax.set_aspect("equal")
ax.text(-0.22, 1.04, "a", transform=ax.transAxes, fontsize=10, fontweight="bold")
ax = axs[1]
for k, lab, c, lw in MODELS:
    s = D[k].values; p, r, _ = precision_recall_curve(y, s); lo, hi = band(pr_on_grid, s)
    ax.fill_between(gx, lo, hi, color=c, alpha=0.18, lw=0); ax.plot(r, p, color=c, lw=lw)
ax.axhline(y.mean(), ls=":", color="#b0b0b0", lw=0.8); ax.set_xlabel("Recall"); ax.set_ylabel("Precision")
ax.set_xlim(-0.02, 1.02); ax.set_ylim(0.3, 1.02); ax.text(-0.22, 1.04, "b", transform=ax.transAxes, fontsize=10, fontweight="bold")
h, l = axs[0].get_legend_handles_labels(); fig.legend(h, l, frameon=False, loc="lower center", ncol=3, bbox_to_anchor=(0.5, -0.07), handlelength=1.8, columnspacing=1.6)
fig.tight_layout(); fig.savefig(FIG/"fig1_roc_pr.png", bbox_inches="tight"); fig.savefig(FIG/"fig1_roc_pr.pdf", bbox_inches="tight"); plt.close(fig)
# ---------------- Figure 2: reliability (equal-frequency bins, Wilson 95% intervals on the observed fraction)
def wilson(k, n, z=1.96):
    p = k/n; d = 1 + z*z/n; c = (p + z*z/(2*n))/d; h = z*np.sqrt(p*(1-p)/n + z*z/(4*n*n))/d; return c-h, c+h
fig, ax = plt.subplots(figsize=(3.3, 3.3))
ax.plot([0, 1], [0, 1], color="#555555", lw=0.8, ls="--", label="Perfect calibration")
for k, lab, c, lw in MODELS[:2]:
    s = D[k].values; order = np.argsort(s); bins = np.array_split(order, 10)
    mp = np.array([s[b].mean() for b in bins]); ob = np.array([y[b].mean() for b in bins]); n = np.array([len(b) for b in bins])
    lo, hi = wilson(ob*n, n); ax.errorbar(mp, ob, yerr=[ob-lo, hi-ob], color=c, lw=lw, marker="o", ms=3.5, capsize=1.5, elinewidth=0.7, label=lab)
ax.set_xlabel("Mean predicted probability of solubility"); ax.set_ylabel("Observed fraction soluble"); ax.set_xlim(-0.03, 1.03); ax.set_ylim(-0.03, 1.03); ax.set_aspect("equal")
ax.legend(frameon=False, loc="upper left", handlelength=1.6)
fig.tight_layout(); fig.savefig(FIG/"fig2_calibration.png", bbox_inches="tight"); fig.savefig(FIG/"fig2_calibration.pdf", bbox_inches="tight"); plt.close(fig)
print("figures written")

