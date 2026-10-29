
"""Naive biophysical baseline (pre-registered): logistic regression on length, GRAVY, net charge at pH 7, charged fraction, aromatic fraction."""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
import numpy as np, pandas as pd
from Bio.SeqUtils.ProtParamData import kd
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
FEATS = ["length", "gravy", "net_charge_pH7", "frac_charged", "frac_aromatic"]
def featurize(seqs):
    rows = []
    for s in seqs:
        pa = ProteinAnalysis(s); L = len(s)
        rows.append(dict(length=L, gravy=pa.gravy(), net_charge_pH7=pa.charge_at_pH(7.0),
                         frac_charged=sum(s.count(c) for c in "DEKR")/L, frac_aromatic=sum(s.count(c) for c in "FWY")/L))
    return pd.DataFrame(rows)[FEATS]
def fit_baseline(train_seqs, train_y):
    m = make_pipeline(StandardScaler(), LogisticRegression(max_iter=1000))   # sklearn default L2, C=1
    return m.fit(featurize(train_seqs), np.asarray(train_y))
def predict_baseline(m, seqs): return m.predict_proba(featurize(seqs))[:, 1]
