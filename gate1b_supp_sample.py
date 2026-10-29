import sys, os, json
sys.path.insert(0, os.path.dirname(__file__))
import numpy as np, pandas as pd, lib
tr = lib.load_uesolds_train(lib.DATA/"plmsol_train_dataset.fasta"); te = lib.load_uesolds(lib.DATA/"plmsol_test_dataset.fasta")
print("train N", len(tr), "label mean", round(tr.label.mean(), 3))
ov = tr.md5.isin(set(te.md5)); tl = dict(zip(te.md5, te.label))
conc = (tr[ov].label.values == tr[ov].md5.map(tl).values)
ev = pd.read_csv(lib.DATA/"eval_sample.csv"); ev_in_train = ev.md5.isin(set(tr.md5))
info = dict(train_n=len(tr), train_rows_with_exact_seq_in_test=int(ov.sum()), of_which_label_concordant=int(conc.sum()),
            test_seqs_with_exact_dup_in_train=int(te.md5.isin(set(tr.md5)).sum()), eval_sample_seqs_with_exact_dup_in_train=int(ev_in_train.sum()))
print(info); json.dump(info, open(lib.RES/"benchmark_train_test_overlap.json", "w"), indent=1)
rng = np.random.default_rng(lib.SEED + 1)
pool = tr[(~ov) & (tr.length >= 50) & (tr.length <= 400)]
s = pd.concat([pool[pool.label == k].sample(50, random_state=int(rng.integers(1e9))) for k in (0, 1)]).reset_index(drop=True)
s.to_csv(lib.DATA/"gate1b_supp_train_sample.csv", index=False); print(s.groupby("label").length.agg(["count", "mean"]))
