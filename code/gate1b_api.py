
import sys, json, time
sys.path.insert(0, __import__("os").path.dirname(__file__))
import numpy as np, pandas as pd
from scipy.stats import mannwhitneyu
import lib
c = pd.read_csv(lib.DATA/"controls.csv")
api = lib.AikiSolAPI(delay=1.0)
rows = []
t0 = time.time()
for i, r in c.iterrows():
    out = api.lookup(r.seq, tag="gate1b")
    out.update(accession=r.accession, name=r["name"], cls=r.cls, label=r.label, length=r.length)
    rows.append(out)
res = pd.DataFrame(rows); res.to_csv(lib.RES/"gate1b_api_predictions.csv", index=False)
print("errors:", int(res.get("error", pd.Series(dtype=object)).notna().sum()), "elapsed s:", round(time.time()-t0))
cols = ["prob_marginal","mean_prob"] + ["prob_"+r for r in lib.REGIMES]
summ = {}
for k in cols:
    a = res.loc[res.label==1, k].astype(float); b = res.loc[res.label==0, k].astype(float)
    summ[k] = dict(mean_soluble=float(a.mean()), mean_insoluble=float(b.mean()), diff=float(a.mean()-b.mean()),
                   auc=lib.auroc(res.label, res[k].astype(float)), mwu_p_one_sided=float(mannwhitneyu(a, b, alternative="greater").pvalue))
passed = bool(summ["prob_marginal"]["mean_soluble"] > summ["prob_marginal"]["mean_insoluble"])
out = dict(n_soluble=int((res.label==1).sum()), n_insoluble=int((res.label==0).sum()), summary=summ, pass_gate_marginal_direction=passed, api_requests_total=lib.api_count())
json.dump(out, open(lib.RES/"gate1b_api.json","w"), indent=1)
print(pd.DataFrame(summ).T.round(4).to_string()); print("PASS (marginal soluble > insoluble):", passed, "| API requests so far:", lib.api_count())
