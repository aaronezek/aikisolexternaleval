
import sys, json, hashlib
sys.path.insert(0, __import__("os").path.dirname(__file__))
import numpy as np, pandas as pd, lib
te = lib.load_uesolds(lib.DATA/"plmsol_test_dataset.fasta")
rng = np.random.default_rng(lib.SEED)               # pre-registered: default_rng(20261028), 250 per class
parts = []
for k in (0, 1):
    idx = te.index[te.label == k].to_numpy()
    parts.append(te.loc[np.sort(rng.choice(idx, 250, replace=False))])
ev = pd.concat(parts).sort_index().reset_index().rename(columns={"index": "test_row"})
ev.to_csv(lib.DATA/"eval_sample.csv", index=False)
print("N", len(ev), ev.label.value_counts().to_dict()); print(ev.length.describe().round(1).to_dict())
print(ev.source.value_counts().to_dict())
print("label rate by source:", ev.groupby("source").label.agg(["mean","count"]).round(2).T.to_dict())
# leakage list format
lines = open(lib.DATA/"zenodo/apache_tier__canonical_147k_md5.txt").read().split("\n")
print("md5 list lines:", len(lines), lines[:3])
