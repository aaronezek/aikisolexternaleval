
"""Local CPU scorer wrapping the UNMODIFIED published aikisol package (source commit d8fcbf86) + Zenodo canonical checkpoint.
Usage: python local_aikisol.py in.csv out.csv [seq_col] [--no-normalize]
Only change vs the package defaults: torch.load(mmap=True) to cut peak RAM on a 16 GB machine (does not change weights/math)."""
import sys, os, time
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/"aikisol_src"))
import numpy as np, pandas as pd, torch
torch.set_num_threads(os.cpu_count() or 8)
_orig_load = torch.load
def _mmap_load(*a, **k):
    try: return _orig_load(*a, mmap=True, **k)
    except Exception as e:
        print("mmap load failed, falling back:", type(e).__name__, flush=True); return _orig_load(*a, **k)
torch.load = _mmap_load
from aikisol import Aikisol
from aikisol.inference import STRINGENCY_LABELS
CKPT = ROOT/"data/zenodo/aikisol_v2_canonical_147k_full.pt"

def load_model():
    return Aikisol.from_pretrained(checkpoint=CKPT, device="cpu")

def score(model, seqs, normalize=True, chunk=32, batch_size=8, log_every=1):
    seqs = list(seqs); order = np.argsort([len(s) for s in seqs], kind="stable")
    cols = {k: np.full(len(seqs), np.nan) for k in ["mean_prob","prob_marginal"] + ["prob_"+r for r in STRINGENCY_LABELS]}
    norm = [None]*len(seqs); t0 = time.time()
    for ci in range(0, len(order), chunk):
        ix = order[ci:ci+chunk]
        out = model.predict([seqs[i] for i in ix], normalize=normalize, batch_size=batch_size)
        cols["mean_prob"][ix] = out.mean_prob; cols["prob_marginal"][ix] = out.prob_marginal
        for r in STRINGENCY_LABELS: cols["prob_"+r][ix] = out.per_stringency[r]
        for j, i in enumerate(ix): norm[i] = out.normalized_sequences[j]
        if (ci//chunk) % log_every == 0: print(f"  scored {min(ci+chunk,len(seqs))}/{len(seqs)}  {time.time()-t0:.0f}s", flush=True)
    d = pd.DataFrame(cols); d["normalized_sequence"] = norm; return d

if __name__ == "__main__":
    inp, outp = sys.argv[1], sys.argv[2]
    col = sys.argv[3] if len(sys.argv) > 3 and not sys.argv[3].startswith("--") else "seq"
    df = pd.read_csv(inp); t0 = time.time(); m = load_model(); print(f"model loaded in {time.time()-t0:.0f}s", flush=True)
    res = score(m, df[col].tolist(), normalize="--no-normalize" not in sys.argv)
    pd.concat([df.reset_index(drop=True), res], axis=1).to_csv(outp, index=False); print("done", time.time()-t0, flush=True)
