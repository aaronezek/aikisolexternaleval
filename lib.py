
"""Shared, verified-once library for the Aiki-Sol external evaluation. Reuse; do not fork."""
import os, re, json, time, hashlib, math
from pathlib import Path
import numpy as np, pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DATA, RES, LOGS = ROOT/"data", ROOT/"results", ROOT/"logs"
for _d in (DATA, RES, LOGS): _d.mkdir(exist_ok=True)
SEED = 20261028
AA20 = "ACDEFGHIKLMNPQRSTVWY"

# ---------------------------------------------------------------- data loading
def read_fasta(path):
    recs, cur = [], None
    for line in open(path, encoding="utf-8"):
        line = line.rstrip("\n")
        if line.startswith(">"):
            if cur: recs.append(cur)
            cur = {"header": line[1:], "seq": ""}
        else:
            cur["seq"] += line.strip()
    if cur: recs.append(cur)
    return pd.DataFrame(recs)

def load_uesolds(path):
    """UESolDS FASTA (PLM_Sol). Header: '<md5> <source-id> A-<label>'; label 1 = soluble, 0 = insoluble."""
    d = read_fasta(path)
    parts = d.header.str.split(" ")
    assert (parts.str.len() == 3).all(), "unexpected header format"
    d["md5"] = parts.str[0]; d["source_id"] = parts.str[1]
    d["label"] = parts.str[2].str.extract(r"^A-(\d)$")[0].astype(int)
    d["source"] = d.source_id.str.extract(r"^([A-Za-z]+)")[0].where(~d.source_id.str.match(r"^[0-9][A-Za-z0-9]{3}_"), "PDB")
    d["length"] = d.seq.str.len()
    assert (d.md5 == d.seq.map(lambda s: hashlib.md5(s.encode()).hexdigest())).all(), "md5 header != md5(sequence)"
    return d.drop(columns="header")

# ---------------------------------------------------------------- SWI (Bhandari et al. 2020; weights copied from reference swi.py)
SWI_W = {'A':0.8356471476582918,'C':0.5208088354857734,'E':0.9876987431418378,'D':0.9079044671339564,'G':0.7997168496420723,
 'F':0.5849790194237692,'I':0.6784124413866582,'H':0.8947913996466419,'K':0.9267104557513497,'M':0.6296623675420369,
 'L':0.6554221515081433,'N':0.8597433107431216,'Q':0.789434648348208,'P':0.8235328714705341,'S':0.7440908318492778,
 'R':0.7712466317693457,'T':0.8096922697856334,'W':0.6374678690957594,'V':0.7357837119163659,'Y':0.6112801822947587}
SWI_A, SWI_B = 81.0581, -62.7775
def swi(seq):
    return float(np.mean([SWI_W[c] for c in seq]))
def swi_prob(seq):
    return 1.0/(1.0+math.exp(-(SWI_A*swi(seq)+SWI_B)))

# ---------------------------------------------------------------- metrics
from sklearn.metrics import (roc_auc_score, average_precision_score, matthews_corrcoef, brier_score_loss, roc_curve)
from sklearn.linear_model import LogisticRegression

def auroc(y, s): return float(roc_auc_score(y, s))
def auprc(y, s): return float(average_precision_score(y, s))

def binary_metrics(y, s, thr=0.5):
    y = np.asarray(y); pred = (np.asarray(s) >= thr).astype(int)
    tp = int(((pred==1)&(y==1)).sum()); tn = int(((pred==0)&(y==0)).sum())
    fp = int(((pred==1)&(y==0)).sum()); fn = int(((pred==0)&(y==1)).sum())
    return dict(thr=float(thr), acc=(tp+tn)/len(y), sens=tp/max(tp+fn,1), spec=tn/max(tn+fp,1),
                mcc=float(matthews_corrcoef(y, pred)) if len(set(pred))>1 else 0.0, tp=tp, tn=tn, fp=fp, fn=fn)

def youden_threshold(y, s):
    fpr, tpr, thr = roc_curve(y, s); j = tpr - fpr
    return float(thr[int(np.argmax(j))])

def _logit(p, eps=1e-6):
    p = np.clip(np.asarray(p, float), eps, 1-eps); return np.log(p/(1-p))

def calibration_stats(y, p, n_bins=10):
    """Brier, ECE (equal-frequency bins), calibration intercept/slope (logistic recalibration on logit(p))."""
    y = np.asarray(y); p = np.asarray(p, float)
    brier = float(brier_score_loss(y, p))
    order = np.argsort(p); bins = np.array_split(order, n_bins)
    ece = float(sum(len(b)/len(y) * abs(y[b].mean() - p[b].mean()) for b in bins))
    lr = LogisticRegression(penalty=None, max_iter=1000).fit(_logit(p).reshape(-1,1), y)
    return dict(brier=brier, ece=ece, cal_intercept=float(lr.intercept_[0]), cal_slope=float(lr.coef_[0][0]),
                mean_pred=float(p.mean()), obs_rate=float(y.mean()))

def reliability_bins(y, p, n_bins=10):
    y = np.asarray(y); p = np.asarray(p, float); order = np.argsort(p); out = []
    for b in np.array_split(order, n_bins):
        out.append((float(p[b].mean()), float(y[b].mean()), len(b)))
    return out

# ---------------------------------------------------------------- bootstrap
def boot_indices(n, B=2000, seed=SEED):
    rng = np.random.default_rng(seed); return rng.integers(0, n, size=(B, n))

def bootstrap_ci(fn, y, *scores, B=2000, seed=SEED, alpha=0.05):
    """Percentile bootstrap of fn(y, *scores). Resamples that contain one class only are skipped (reported in n_valid)."""
    y = np.asarray(y); scores = [np.asarray(s) for s in scores]
    idx = boot_indices(len(y), B, seed); vals = []
    for ix in idx:
        yy = y[ix]
        if yy.min() == yy.max(): continue
        vals.append(fn(yy, *[s[ix] for s in scores]))
    vals = np.array(vals)
    return dict(point=float(fn(y, *scores)), lo=float(np.quantile(vals, alpha/2)), hi=float(np.quantile(vals, 1-alpha/2)), n_valid=int(len(vals)))

# ---------------------------------------------------------------- Aiki-Sol REST wrapper (every HTTP call is counted and logged)
API_BASE = "https://aikium--aikisol-landing-page.modal.run"
API_LOG = LOGS/"api_requests.jsonl"
API_CAP = 1500
REGIMES = ["3000g_10min", "6000g", "32000g", "eSol_21600g_30min", "100000g"]

def api_count():
    return sum(1 for _ in open(API_LOG, encoding="utf-8")) if API_LOG.exists() else 0

class AikiSolAPI:
    def __init__(self, delay=1.0, cap=API_CAP, timeout=120):
        import requests
        self.s = requests.Session(); self.delay = delay; self.cap = cap; self.timeout = timeout; self._last = 0.0
    def _call(self, method, path, tag="", **kw):
        if api_count() >= self.cap: raise RuntimeError(f"API request cap {self.cap} reached")
        wait = self.delay - (time.time() - self._last)
        if wait > 0: time.sleep(wait)
        t0 = time.time(); status = None; text = ""; err = None
        try:
            r = self.s.request(method, API_BASE + path, timeout=self.timeout, **kw)
            status, text = r.status_code, r.text
        except Exception as e:
            err = f"{type(e).__name__}: {str(e)[:200]}"
        self._last = time.time()
        rec = dict(ts=time.strftime("%Y-%m-%dT%H:%M:%S"), method=method, path=path, tag=tag,
                   req=(kw.get("params") or kw.get("json")), status=status, latency_s=round(self._last - t0, 3),
                   response=text[:20000], error=err)
        with open(API_LOG, "a", encoding="utf-8") as fh: fh.write(json.dumps(rec) + "\n")
        return status, text, err
    def lookup_raw(self, seq, tag=""):
        return self._call("GET", "/lookup_protein", tag=tag, params={"seq": seq})
    def predict_fasta_raw(self, fasta, fmt="json", tag=""):
        return self._call("POST", "/predict_fasta", tag=tag, json={"fasta": fasta, "format": fmt})

    def lookup(self, seq, tag="", retries=2):
        """Score one sequence via GET /lookup_protein. Returns flat dict (mean_prob, prob_marginal, prob_<regime>..., normalized_sequence)
        or {'error': ...}. Retries (each counted) only on throttling / server / network errors."""
        st = err = tx = None
        for k in range(retries + 1):
            st, tx, err = self.lookup_raw(seq, tag=tag)
            if st == 200:
                js = json.loads(tx)
                out = dict(mean_prob=js["mean_prob"], prob_marginal=js["prob_marginal"], normalized_sequence=js.get("normalized_sequence"))
                for r in REGIMES: out["prob_" + r] = js["per_stringency"][r]
                return out
            if st in (429, 500, 502, 503, 504) or err:
                time.sleep(5 * (k + 1)); continue
            return dict(error=f"HTTP {st}: {(tx or '')[:300]}")
        return dict(error=f"HTTP {st} {err} {(tx or '')[:200]}")


def load_uesolds_train(path):
    """UESolDS TRAIN FASTA. Header: '<id> A-<label>' (2 fields; no md5/source). label 1 = soluble. Added for the supplementary Gate-1b check;
    the verified test-split loader load_uesolds() is unchanged."""
    d = read_fasta(path)
    parts = d.header.str.split(" ")
    assert (parts.str.len() == 2).all(), "unexpected train header format"
    d["source_id"] = parts.str[0]
    d["label"] = parts.str[1].str.extract(r"^A-(\d)$")[0].astype(int)
    d["md5"] = d.seq.map(lambda s: hashlib.md5(s.encode()).hexdigest())
    d["length"] = d.seq.str.len()
    return d.drop(columns="header")
