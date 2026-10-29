
"""Milestones 3-4: baselines + evaluation. PRIMARY = AUROC of Aiki-Sol prob_marginal on the 500-sequence UESolDS sample (pre-registered).
Everything else is secondary/exploratory and labelled as such in the output tables."""
import sys, os, json
sys.path.insert(0, os.path.dirname(__file__))
import numpy as np, pandas as pd, lib, baselines
RES = lib.RES
ev = pd.read_csv(lib.DATA/"eval_sample.csv")
pr = pd.read_csv(RES/"eval_local_predictions.csv")
assert len(pr) == len(ev) and (pr.md5 == ev.md5).all()
leak = pd.read_csv(RES/"leakage_exact_eval.csv"); assert (leak.md5 == ev.md5).all()
ev["exact_in_train"] = leak.exact_in_train.values
y = ev.label.values
OUT = {"n": len(ev), "n_soluble": int(y.sum())}
SC = {"aikisol_marginal": pr.prob_marginal.values, "aikisol_mean5": pr.mean_prob.values}
for r in lib.REGIMES: SC["aikisol_"+r] = pr["prob_"+r].values

# ---- baselines
ev["swi_prob"] = ev.seq.map(lib.swi_prob); SC["swi"] = ev.swi_prob.values
tr_path = lib.DATA/"plmsol_train_dataset.fasta"
if tr_path.exists():
    tr = lib.load_uesolds_train(tr_path); assert not set(tr.md5) & set(ev.md5); n_all = len(tr); tr = tr[tr.seq.str.fullmatch("[ACDEFGHIKLMNPQRSTVWY]+")]; OUT["train_noncanonical_dropped"] = n_all - len(tr)
    lr = baselines.fit_baseline(tr.seq.tolist(), tr.label.values); OUT["lr_training_data"] = f"UESolDS train split (N={len(tr)}), as pre-registered"
    SC["logreg5"] = baselines.predict_baseline(lr, ev.seq.tolist())
else:   # documented deviation: fit on the 3,500 test-split sequences NOT in the evaluation sample
    te = lib.load_uesolds(lib.DATA/"plmsol_test_dataset.fasta"); rest = te[~te.md5.isin(ev.md5)]
    lr = baselines.fit_baseline(rest.seq.tolist(), rest.label.values); OUT["lr_training_data"] = f"DEVIATION: UESolDS train FASTA unreachable (network); fit on the {len(rest)} test-split sequences outside the evaluation sample"
    SC["logreg5"] = baselines.predict_baseline(lr, ev.seq.tolist())
coef = dict(zip(baselines.FEATS, lr[-1].coef_[0].round(3))); OUT["lr_coefficients_standardised"] = {k: float(v) for k, v in coef.items()}
SC["majority"] = np.full(len(y), y.mean())

def row(name, s, thr=0.5):
    a = lib.bootstrap_ci(lib.auroc, y, s); p = lib.bootstrap_ci(lib.auprc, y, s)
    d = dict(model=name, auroc=a["point"], auroc_lo=a["lo"], auroc_hi=a["hi"], auprc=p["point"], auprc_lo=p["lo"], auprc_hi=p["hi"])
    for tag, t in (("t0.5", thr), ("tYouden", lib.youden_threshold(y, s))):
        bm = lib.binary_metrics(y, s, t); d[f"thr_{tag}"] = t
        for m in ("mcc", "acc", "sens", "spec"):
            f = (lambda m, t: (lambda yy, ss: lib.binary_metrics(yy, ss, t)[m]))(m, t)
            c = lib.bootstrap_ci(f, y, s); d[f"{m}_{tag}"] = c["point"]; d[f"{m}_{tag}_lo"] = c["lo"]; d[f"{m}_{tag}_hi"] = c["hi"]
    return d
rows = [row(k, v) for k, v in SC.items() if k != "majority"]
maj = dict(model="majority", auroc=0.5, auroc_lo=0.5, auroc_hi=0.5, auprc=float(y.mean()), acc_t0p5=float(max(y.mean(), 1-y.mean())))
tab = pd.DataFrame(rows + [maj]); tab.to_csv(RES/"table_all_models.csv", index=False)
OUT["primary"] = dict(metric="AUROC prob_marginal", **lib.bootstrap_ci(lib.auroc, y, SC["aikisol_marginal"]))

# ---- paired differences (AUROC) vs baselines, same bootstrap indices
diffs = {}
for b in ("logreg5", "swi"):
    for a in ("aikisol_marginal",) + tuple("aikisol_"+r for r in lib.REGIMES):
        diffs[f"{a} - {b}"] = lib.bootstrap_ci(lambda yy, s1, s2: lib.auroc(yy, s1) - lib.auroc(yy, s2), y, SC[a], SC[b])
OUT["paired_auroc_diffs"] = diffs

# ---- calibration (all regime outputs + marginal; primary = marginal)
cal = {k: lib.calibration_stats(y, SC[k]) for k in SC if k.startswith("aikisol_") or k in ("swi", "logreg5")}
for k in ("aikisol_marginal", "logreg5"):
    for m in ("brier", "ece"):
        cal[k][m+"_ci"] = lib.bootstrap_ci(lambda yy, pp, m=m: lib.calibration_stats(yy, pp)[m], y, SC[k], B=1000)
OUT["calibration"] = cal; OUT["reliability_marginal"] = lib.reliability_bins(y, SC["aikisol_marginal"])
OUT["reliability_logreg5"] = lib.reliability_bins(y, SC["logreg5"])
# prior-shift note: the benchmark is 50% soluble; mean predicted vs observed
OUT["mean_pred_marginal"] = float(SC["aikisol_marginal"].mean())

# ---- length (secondary)
ev["len_tercile"] = pd.qcut(ev.length, 3, labels=["short", "medium", "long"]); L = {}
for g, d in ev.groupby("len_tercile", observed=True):
    ix = d.index.values; L[g] = dict(n=len(ix), len_range=[int(d.length.min()), int(d.length.max())],
        **{m: lib.bootstrap_ci(lib.auroc, y[ix], SC[m][ix]) for m in ("aikisol_marginal", "logreg5", "swi")})
OUT["by_length_tercile"] = L
from scipy.stats import spearmanr
OUT["spearman_score_vs_length"] = float(spearmanr(ev.length, SC["aikisol_marginal"])[0])
# ---- leakage strata (secondary, central to interpretation)
S = {}
for name, mask in (("exact_match_to_training", ev.exact_in_train.values), ("no_exact_match", ~ev.exact_in_train.values)):
    ix = np.where(mask)[0]; S[name] = dict(n=len(ix), n_soluble=int(y[ix].sum()),
        **{m: lib.bootstrap_ci(lib.auroc, y[ix], SC[m][ix]) for m in ("aikisol_marginal", "logreg5", "swi")})
OUT["by_leakage_stratum"] = S
# ---- source confound (exploratory): AUROC within source; source-only label-rate predictor (in-sample, an upper bound on what source alone gives)
rate = ev.groupby("source").label.transform("mean").values; OUT["source_only_auroc_in_sample"] = lib.auroc(y, rate)
SS = {}
for src, d in ev.groupby("source"):
    ix = d.index.values
    if min(y[ix].sum(), len(ix)-y[ix].sum()) >= 10: SS[src] = dict(n=len(ix), n_soluble=int(y[ix].sum()), **{m: lib.bootstrap_ci(lib.auroc, y[ix], SC[m][ix]) for m in ("aikisol_marginal", "logreg5")})
OUT["within_source"] = SS
OUT["source_counts"] = ev.source.value_counts().to_dict()
json.dump(OUT, open(RES/"eval_results.json", "w"), indent=1, default=float)
pd.DataFrame({"md5": ev.md5, "label": y, "source": ev.source, "length": ev.length, "exact_in_train": ev.exact_in_train, **{k: v for k, v in SC.items()}}).to_csv(RES/"eval_scores_all_models.csv", index=False)
p = OUT["primary"]; print(f"PRIMARY  AUROC prob_marginal = {p['point']:.3f}  95% CI [{p['lo']:.3f}, {p['hi']:.3f}]  (N={len(y)})")
print(tab[["model","auroc","auroc_lo","auroc_hi","auprc","mcc_t0.5","acc_t0.5"]].round(3).to_string(index=False) if "mcc_t0.5" in tab else tab.iloc[:, :8].round(3).to_string(index=False))


