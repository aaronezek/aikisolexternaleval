
import sys, json, hashlib
sys.path.insert(0, __import__("os").path.dirname(__file__))
import numpy as np, pandas as pd, lib
md5 = set(x.strip() for x in open(lib.DATA/"zenodo/apache_tier__canonical_147k_md5.txt") if x.strip())
print("training-pool md5 entries:", len(md5))
H = lambda s: hashlib.md5(s.encode()).hexdigest()
def variants(s):
    v = {s, "M"+s}
    if s.startswith("M"): v.add(s[1:])
    return v
te = lib.load_uesolds(lib.DATA/"plmsol_test_dataset.fasta"); ev = pd.read_csv(lib.DATA/"eval_sample.csv")
def flag(df): return df.seq.map(lambda s: any(H(v) in md5 for v in variants(s)))
te["exact_in_train"] = flag(te); ev["exact_in_train"] = flag(ev)
print("full test split (4000): exact matches to Aiki-Sol training pool:", int(te.exact_in_train.sum()), f"({te.exact_in_train.mean():.2%})")
print("eval sample (500):", int(ev.exact_in_train.sum()), f"({ev.exact_in_train.mean():.2%})")
print(ev.groupby("source").exact_in_train.agg(["sum","count"]).T.to_dict())
print("by label (eval):", ev.groupby("label").exact_in_train.sum().to_dict())
ev[["md5","exact_in_train"]].to_csv(lib.RES/"leakage_exact_eval.csv", index=False)
json.dump(dict(test_4000_exact=int(te.exact_in_train.sum()), eval_500_exact=int(ev.exact_in_train.sum()),
               md5_list_n=len(md5), note="exact md5 of raw, M-prefixed, and M-stripped sequence vs canonical-147K md5 list"), open(lib.RES/"leakage_exact.json","w"), indent=1)
