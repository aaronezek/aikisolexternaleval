
import sys, json, os
sys.path.insert(0, os.path.dirname(__file__))
import numpy as np, pandas as pd, lib
a = pd.read_csv(lib.RES/"gate1b_api_predictions.csv"); y = a.label.values; out = {}
for k in ["prob_marginal","mean_prob"] + ["prob_"+r for r in lib.REGIMES]:
    out[k] = lib.bootstrap_ci(lib.auroc, y, a[k].values)
print(pd.DataFrame(out).T[["point","lo","hi"]].round(3).to_string())
json.dump(out, open(lib.RES/"gate1b_control_auc_ci.json","w"), indent=1)
print("train fasta present:", (lib.DATA/"plmsol_train_dataset.fasta").exists(), "| API requests so far:", lib.api_count())
