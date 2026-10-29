
import sys, json
sys.path.insert(0, __import__("os").path.dirname(__file__))
import numpy as np, pandas as pd
import lib
d = lib.load_uesolds(lib.DATA/"plmsol_test_dataset.fasta")
n, npos = len(d), int(d.label.sum())
print("N =", n, "soluble =", npos, "insoluble =", n-npos)
assert n == 4000 and npos == 2000, "split size/balance mismatch vs published"
assert d.seq.str.fullmatch("[ACDEFGHIKLMNPQRSTVWY]+").all()
d["swi"] = d.seq.map(lib.swi); d["swi_prob"] = d.seq.map(lib.swi_prob)
y = d.label.values
auc = lib.bootstrap_ci(lib.auroc, y, d.swi_prob.values)
bm = lib.binary_metrics(y, d.swi_prob.values, 0.5)
pub = dict(auc=0.5597, acc=0.5423, mcc=0.0971, sens=0.7905, spec=0.2938)
res = dict(n=n, n_soluble=npos, auc=auc, metrics_at_0p5=bm, published=pub,
           abs_diff_auc=abs(auc["point"]-pub["auc"]), tolerance=0.05,
           pass_gate=bool(abs(auc["point"]-pub["auc"]) <= 0.05 and n == 4000 and npos == 2000))
print(json.dumps(res, indent=1))
d[["md5","source_id","source","label","length","swi","swi_prob"]].to_csv(lib.RES/"gate1a_swi_predictions.csv", index=False)
json.dump(res, open(lib.RES/"gate1a.json","w"), indent=1)
