import sys, os, json
sys.path.insert(0, os.path.dirname(__file__))
import numpy as np, pandas as pd, lib
from scipy.stats import mannwhitneyu
d = pd.read_csv(lib.RES/"gate1b_supp_local_predictions.csv"); y = d.label.values
print("N", len(d), "soluble", int(y.sum()))
out = {}
for k in ["prob_marginal", "mean_prob"] + ["prob_"+r for r in lib.REGIMES]:
    c = lib.bootstrap_ci(lib.auroc, y, d[k].values, B=1000)
    a, b = d.loc[y==1, k], d.loc[y==0, k]
    out[k] = dict(auc=c["point"], lo=c["lo"], hi=c["hi"], mean_sol=float(a.mean()), mean_insol=float(b.mean()), mwu_p=float(mannwhitneyu(a, b, alternative="greater").pvalue))
t = pd.DataFrame(out).T.round(3); print(t.to_string())
json.dump(out, open(lib.RES/"gate1b_supp.json", "w"), indent=1)
