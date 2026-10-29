"""M5 probes against the HOSTED API (exploratory). Reuses lib.AikiSolAPI (every request counted/logged)."""
import sys, os, json, time
sys.path.insert(0, os.path.dirname(__file__))
import numpy as np, pandas as pd, lib
T = pd.read_csv(lib.RES/"probe_tags_local.csv")
api = lib.AikiSolAPI(delay=1.0)
bases = T[T.variant == "base"].drop_duplicates("md5")
# ---- 1. tag probe on 6 sequences x 7 variants (42 requests)
sub = bases.md5.tolist()[:6]; rows = []
for _, r in T[T.md5.isin(sub)].iterrows():
    o = api.lookup(r.seq_in, tag="probe_tags"); o.update(md5=r.md5, variant=r.variant, label=r.label, seq_in_len=len(r.seq_in)); rows.append(o)
A = pd.DataFrame(rows); A.to_csv(lib.RES/"probe_tags_api.csv", index=False)
# ---- 2. determinism: 10 sequences x 2 identical calls (20 requests)
d10 = bases.seq_in.tolist()[:10]; det = []
for rep in (1, 2):
    for s in d10:
        o = api.lookup(s, tag=f"probe_determinism_rep{rep}"); o["rep"] = rep; o["seq_len"] = len(s); det.append(o)
D = pd.DataFrame(det); D.to_csv(lib.RES/"probe_determinism_api.csv", index=False)
cols = ["prob_marginal", "mean_prob"] + ["prob_" + r for r in lib.REGIMES]
a = D[D.rep == 1].reset_index(drop=True)[cols].astype(float); b = D[D.rep == 2].reset_index(drop=True)[cols].astype(float)
det_res = dict(n=len(d10), max_abs_diff=float((a - b).abs().values.max()), identical=bool((a == b).all().all()))
print("determinism (API):", det_res)
# ---- 3. robustness (16 requests); raw HTTP status + body excerpt recorded
b0 = bases.seq_in.iloc[0]; b1 = bases.seq_in.iloc[1]; rng = np.random.default_rng(lib.SEED + 7)
rob_in = {"empty": "", "single_M": "M", "dipeptide": "MK", "5mer": "MKLIT", "lowercase": b0.lower(), "with_X": b0[:60] + "XXXX" + b0[60:], "with_B": b0[:60] + "B" + b0[61:],
          "with_Z": b0[:60] + "Z" + b0[61:], "with_U": b0[:60] + "U" + b0[61:], "digits_and_star": "MKLIT123VLV*AAK", "whitespace_inside": "MKLITV LVLALLAV", "all_A_300": "A" * 300,
          "poly_H_100": "M" + "H" * 100, "random_300": "".join(rng.choice(list(lib.AA20), 300)), "long_1500_real": (b1 * 8)[:1500], "long_3000_real": (b1 * 16)[:3000]}
rob = []
for k, s in rob_in.items():
    st, tx, err = api.lookup_raw(s, tag="probe_robustness")
    rec = dict(case=k, len_in=len(s), http=st, err=err)
    if st == 200:
        js = json.loads(tx); rec.update(prob_marginal=js["prob_marginal"], norm_len=len(js.get("normalized_sequence") or ""))
    else: rec["body"] = (tx or "")[:200]
    rob.append(rec)
R = pd.DataFrame(rob); R.to_csv(lib.RES/"probe_robustness_api.csv", index=False); print(R.to_string())
json.dump(dict(determinism_api=det_res, api_requests_total=lib.api_count()), open(lib.RES/"probe_api_meta.json", "w"), indent=1)
print("API requests total:", lib.api_count())
