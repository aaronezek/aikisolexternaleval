import sys, os, json
sys.path.insert(0, os.path.dirname(__file__))
import numpy as np, pandas as pd, lib
T = pd.read_csv(lib.RES/"probe_tags_local.csv"); A = pd.read_csv(lib.RES/"probe_tags_api.csv")
PAIRS = [("hisN_pET28", "ctrlN_rnd6"), ("hisC", "ctrlC_rnd6"), ("sumoN", "ctrlN_sumo_scrambled")]
def shifts(df, col="prob_marginal"):
    w = df.pivot(index="md5", columns="variant", values=col); return w.sub(w["base"], axis=0).drop(columns="base")
def summarize(df, name):
    sh = shifts(df); out = {}
    for real, ctrl in PAIRS:
        r, c = sh[real].values, sh[ctrl].values
        f = lambda rr, cc: float(np.mean(np.abs(rr)) - np.mean(np.abs(cc)))
        rng = np.random.default_rng(lib.SEED); bs = [f(r[i], c[i]) for i in (rng.integers(0, len(r), len(r)) for _ in range(2000))]
        out[f"{real} vs {ctrl}"] = dict(n=len(r), mean_shift_real=float(r.mean()), mean_abs_shift_real=float(np.abs(r).mean()), mean_shift_ctrl=float(c.mean()), mean_abs_shift_ctrl=float(np.abs(c).mean()),
                                        diff_abs_real_minus_ctrl=f(r, c), ci=[float(np.quantile(bs, .025)), float(np.quantile(bs, .975))])
    print(f"== {name}"); print(pd.DataFrame(out).T.round(4).to_string()); return out
res = {"local": summarize(T, "local (24 sequences)"), "api": summarize(A, "hosted API (6 sequences)")}
# string-level recovery / retention (local normalizer == package normalize_sequence)
g = T.groupby("variant").agg(n=("md5", "size"), tag_flagged=("had_tag", "mean"), norm_equals_base=("norm_equals_base_norm", "mean"), median_extra_residues=("extra_residues_after_norm", "median"), retained_his6=("retained_his6_after_norm", "mean"))
print(g.round(2).to_string()); res["string_level"] = g.round(3).reset_index().to_dict(orient="records")
# local vs API agreement on shared rows
m = T.merge(A, on=["md5", "variant"], suffixes=("_l", "_a")); res["local_vs_api_max_abs_diff_prob_marginal"] = float((m.prob_marginal_l - m.prob_marginal_a).abs().max()); print("local vs API max |diff|:", res["local_vs_api_max_abs_diff_prob_marginal"], "rows", len(m))
json.dump(res, open(lib.RES/"probe_tags_summary.json", "w"), indent=1)
