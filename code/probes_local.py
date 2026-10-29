
"""M5 probes, local route (exploratory). Tag stripping with length-matched controls, determinism, robustness."""
import sys, os, json
sys.path.insert(0, os.path.dirname(__file__))
import numpy as np, pandas as pd, lib
import local_aikisol as L
from aikisol.sequence_normalization import normalize_sequence
rng = np.random.default_rng(lib.SEED + 5)
ev = pd.read_csv(lib.DATA/"eval_sample.csv"); ctr = pd.read_csv(lib.DATA/"controls.csv")
SUMO = ctr.loc[ctr.accession == "Q12306", "seq"].iloc[0]                       # S. cerevisiae Smt3 (UniProt Q12306)
sumo_scr = "".join(rng.permutation(list(SUMO)))                                   # composition-matched scramble
nonH = [a for a in lib.AA20 if a != "H"]
rnd6 = "".join(rng.choice(nonH, 6)); assert "H" not in rnd6
LEAD_PRE, LEAD_POST = "MGSS", "SSGLVPRGSH"                                        # pET-28a-style N-terminal leader
# 24 evaluation-sample sequences, 80-250 aa, 12 soluble / 12 insoluble (probe sample; not used for any metric)
pool = ev[(ev.length >= 80) & (ev.length <= 250)]
sel = pd.concat([pool[pool.label == k].sample(12, random_state=int(rng.integers(1e9))) for k in (0, 1)]).reset_index(drop=True)
rows = []
for _, r in sel.iterrows():
    b = r.seq; core = b[1:] if b.startswith("M") else b
    V = {"base": b,
         "hisN_pET28":  LEAD_PRE + "HHHHHH" + LEAD_POST + core,        # real His6, N-terminal, with thrombin-site leader
         "ctrlN_rnd6":  LEAD_PRE + rnd6 + LEAD_POST + core,            # same leader, His6 replaced by a random non-His 6-mer
         "hisC":        b + "LEHHHHHH",
         "ctrlC_rnd6":  b + "LE" + rnd6,
         "sumoN":       SUMO + core,                                    # SUMO is NOT a documented stripped tag
         "ctrlN_sumo_scrambled": sumo_scr + core}
    for k, v in V.items():
        n, meta = normalize_sequence(v)
        rows.append(dict(md5=r.md5, label=r.label, variant=k, seq_in=v, norm=n, had_tag=meta.had_tag, tag_type=meta.tag_type))
P = pd.DataFrame(rows)
# string-level recovery: does normalisation turn the tagged construct back into the base normalised sequence?
base_norm = P[P.variant == "base"].set_index("md5").norm
P["norm_equals_base_norm"] = P.apply(lambda r: r.norm == base_norm[r.md5], axis=1)
P["extra_residues_after_norm"] = P.apply(lambda r: len(r.norm) - len(base_norm[r.md5]), axis=1)
print(P.groupby("variant").agg(recovered=("norm_equals_base_norm", "mean"), extra=("extra_residues_after_norm", "median"), flagged_tag=("had_tag", "mean")).round(2))
P["retained_his6_after_norm"] = P.norm.str.contains("HHHHHH"); m = L.load_model(); S = L.score(m, P.norm.tolist(), normalize=False, chunk=24)   # per-sequence normalize_sequence (identical to package), no batch preflight
P = pd.concat([P.reset_index(drop=True), S.drop(columns="normalized_sequence")], axis=1)
P.to_csv(lib.RES/"probe_tags_local.csv", index=False)
# determinism (local): same 10 sequences twice in separate calls, plus alone-vs-batched
d10 = sel.seq.tolist()[:10]; a = L.score(m, d10, chunk=10); b = L.score(m, d10, chunk=10)
c = pd.concat([L.score(m, [s], chunk=1) for s in d10], ignore_index=True)
cols = ["prob_marginal", "mean_prob"] + ["prob_" + r for r in lib.REGIMES]
det = dict(repeat_max_abs_diff=float((a[cols] - b[cols]).abs().values.max()), alone_vs_batched_max_abs_diff=float((a[cols] - c[cols]).abs().values.max()))
print("determinism (local):", det)
# robustness (local)
rob_in = {"empty": "", "single_M": "M", "dipeptide": "MK", "5mer": "MKLIT", "lowercase": sel.seq.iloc[0].lower(), "with_X": sel.seq.iloc[0][:60] + "XXXX" + sel.seq.iloc[0][60:],
          "with_B": sel.seq.iloc[0][:60] + "B" + sel.seq.iloc[0][61:], "with_Z": sel.seq.iloc[0][:60] + "Z" + sel.seq.iloc[0][61:], "with_U": sel.seq.iloc[0][:60] + "U" + sel.seq.iloc[0][61:],
          "digits_and_star": "MKLIT123VLV*AAK", "whitespace_inside": "MKLITV LVLALLAV", "all_A_300": "A" * 300, "poly_H_100": "M" + "H" * 100, "random_300": "".join(rng.choice(list(lib.AA20), 300)),
          "long_1500_real": (sel.seq.iloc[1] * 8)[:1500], "long_3000_real": (sel.seq.iloc[1] * 16)[:3000]}
rob = []
for k, s in rob_in.items():
    try:
        o = L.score(m, [s], chunk=1).iloc[0]; rob.append(dict(case=k, len_in=len(s), ok=True, norm_len=len(o.normalized_sequence), prob_marginal=float(o.prob_marginal), error=None))
    except Exception as e: rob.append(dict(case=k, len_in=len(s), ok=False, norm_len=None, prob_marginal=None, error=f"{type(e).__name__}: {str(e)[:120]}"))
R = pd.DataFrame(rob); print(R.to_string()); R.to_csv(lib.RES/"probe_robustness_local.csv", index=False)
json.dump(dict(determinism_local=det, rnd6=rnd6, sumo_seq=SUMO, sumo_scrambled=sumo_scr, n_probe_sequences=len(sel)), open(lib.RES/"probe_local_meta.json", "w"), indent=1)

