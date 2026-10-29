
import sys, os, json
sys.path.insert(0, os.path.dirname(__file__))
import numpy as np, pandas as pd, lib, baselines
ev = pd.read_csv(lib.DATA/"eval_sample.csv"); pr = pd.read_csv(lib.RES/"eval_local_predictions.csv"); lk = pd.read_csv(lib.RES/"leakage_exact_eval.csv")
assert (pr.md5 == ev.md5).all() and (lk.md5 == ev.md5).all()
y = ev.label.values; ex = lk.exact_in_train.values
te = lib.load_uesolds(lib.DATA/"plmsol_test_dataset.fasta"); rest = te[~te.md5.isin(ev.md5)]
lr = baselines.fit_baseline(rest.seq.tolist(), rest.label.values); lrp = baselines.predict_baseline(lr, ev.seq.tolist())
swi = ev.seq.map(lib.swi_prob).values
S = {"AikiSol prob_marginal (PRIMARY)": pr.prob_marginal.values, "AikiSol mean_prob": pr.mean_prob.values}
for r in lib.REGIMES: S["AikiSol "+r] = pr["prob_"+r].values
S["LogReg-5feat (baseline)"] = lrp; S["SWI (published)"] = swi
out = {}
for k, s in S.items(): out[k] = lib.bootstrap_ci(lib.auroc, y, s)
print("N", len(y), "soluble", int(y.sum()))
for k, v in out.items(): print(f"{k:38s} AUROC {v['point']:.3f} [{v['lo']:.3f}, {v['hi']:.3f}]")
d = {}
for b, sb in (("LogReg", lrp), ("SWI", swi)):
    d["marginal - "+b] = lib.bootstrap_ci(lambda yy, a, c: lib.auroc(yy, a) - lib.auroc(yy, c), y, S["AikiSol prob_marginal (PRIMARY)"], sb)
for k, v in d.items(): print(f"paired diff {k:22s} {v['point']:+.3f} [{v['lo']:+.3f}, {v['hi']:+.3f}]")
for nm, m in (("exact-match to training (n=%d)" % ex.sum(), ex), ("no exact match (n=%d)" % (~ex).sum(), ~ex)):
    ix = np.where(m)[0]; a = lib.bootstrap_ci(lib.auroc, y[ix], S["AikiSol prob_marginal (PRIMARY)"][ix]); b = lib.bootstrap_ci(lib.auroc, y[ix], lrp[ix])
    print(f"stratum {nm:36s} soluble={int(y[ix].sum())}  AikiSol {a['point']:.3f} [{a['lo']:.3f},{a['hi']:.3f}]  LogReg {b['point']:.3f} [{b['lo']:.3f},{b['hi']:.3f}]")
p = S["AikiSol prob_marginal (PRIMARY)"]
print("calibration (marginal):", {k: round(v, 3) for k, v in lib.calibration_stats(y, p).items()})
print("calibration (LogReg):  ", {k: round(v, 3) for k, v in lib.calibration_stats(y, lrp).items()})
print("marginal @0.5:", {k: round(v, 3) for k, v in lib.binary_metrics(y, p, 0.5).items() if k in ("acc","sens","spec","mcc")})
json.dump(dict(auroc=out, paired=d), open(lib.RES/"primary_quick.json", "w"), indent=1)
