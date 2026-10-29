# PRE-SPECIFICATION (written before any Aiki-Sol prediction was obtained)

Written: 2026-10-28T11:42:36-07:00
State of knowledge at writing: I had read the Aiki-Sol landing page, GitHub README/source (commit d8fcbf86), the Zenodo
data-source manifest, and the published PLM_Sol paper (Table 2). I had NOT obtained any Aiki-Sol score for any sequence.
(Exception: none. Sample API calls for format documentation in Milestone 0 are made AFTER this file is hashed.)

## Primary endpoint (one number, fixed in advance)
AUROC of Aiki-Sol `prob_marginal` (the protocol-agnostic output; see 2b reasoning) on the evaluation set below,
with a 95% percentile bootstrap CI (B = 2000 resamples, seed 20261028).
Everything else (other metrics, the five per-regime outputs, mean_prob, baselines' CIs beyond the paired comparison,
leakage-filtered subset, per-source / per-length breakdowns, all probes) is SECONDARY / EXPLORATORY and is labelled so.

## Gate 1a (my pipeline: loader + metric code), using a published predictor
Benchmark: UESolDS independent test set (Zhang et al. 2024, PLM_Sol, Brief. Bioinform. 25(5):bbae404), as distributed in
github.com/Violet969/PLM_Sol embedding_dataset/test_dataset.fasta. Published: 4,000 sequences, 2,000 soluble / 2,000 insoluble.
Published predictor: SWI (Bhandari et al. 2020), reference code Gardner-BinfLab/SoDoPE_paper_2020 SWI/swi.py (commit c28f26d9).
Published value on this split (PLM_Sol Table 2, T=0.5): AUC 0.5597, accuracy 0.5423, MCC 0.0971, sensitivity 0.7905, specificity 0.2938.
PASS: loaded N = 4000 with 2000/2000 balance AND |my AUC - 0.5597| <= 0.05.
If FAIL, diagnose in this order: label polarity, wrong split, metric definition, preprocessing.
SoluProt (standalone needs TMHMM + USEARCH, licensed Linux binaries) and NetSolP (5.6 GB bundle) were judged infeasible on this machine
and are NOT re-run; their published numbers may be quoted only as 'cited from PLM_Sol Table 2, full 4,000-sequence split'.

## Gate 1b (my reading of the Aiki-Sol API / output direction)
~20 well-documented soluble-in-E.-coli proteins vs ~20 documented aggregation-prone / inclusion-body proteins (list frozen in code/controls.py
before calling the API). PASS: mean `prob_marginal` of the soluble set > mean of the insoluble set (direction implied by the docs:
output = P(soluble)). If reversed -> suspect my output reading; if indistinguishable -> debug the call, do not report a null.
Also checked (secondary): all five per-regime outputs and API-vs-local-package agreement.

## Evaluation dataset (2a)
UESolDS independent test set (same source as 1a; the 2024 benchmark whose test split was built to share <=25% identity with the
training sets of seven earlier predictors). It is an ESTABLISHED benchmark -> leakage into Aiki-Sol training cannot be ruled out
(Aiki-Sol's pool contains eSol, NESG, PSI:Biology, MCSG/CSGID/NYSGRC TargetTrack, DeepSol, DeepSoluE, ProgSol). Therefore the result is
reported as an upper bound on true generalisation, stated in the abstract. Leakage is audited (exact MD5 match + k-mer/alignment identity
to the public canonical-147K training CSV) and reported; a leakage-filtered subset is exploratory.
Sample: because CPU-only inference of ESM-2 650M is ~seconds per sequence, N = 500, stratified random sample of 250 soluble + 250
insoluble from the 4,000 (numpy default_rng seed 20261028), drawn before scoring. Label 1 = soluble.
Gate-1a and the evaluation set are the SAME benchmark (different sizes: 4,000 for 1a, 500-subsample for evaluation) -- stated throughout.

## Protocol match (2b) -- decided before scoring
Benchmark labels come from heterogeneous structural-genomics pipelines (JCSG, MCSG, SGPP, SECSG, NYSGRC, NESG, SSGCID, ...) plus
PDB-deposition-derived 'soluble' calls; no single centrifugation speed or cell-free condition defines them. None of the five regimes
(3K/6K/32K/100K x g in vivo, PURE 21,600 x g cell-free) plausibly matches the whole set. Per Aiki-Sol's documentation the dedicated
marginal output (6th head) is 'the recommended single number when the downstream protocol is unspecified'. PRIMARY OUTPUT = `prob_marginal`.
All five regimes + mean_prob reported in a secondary table; any regime that beats the primary is an observation about protocol sensitivity,
not the headline. A protocol mismatch is MY limitation, not a flaw of the tool, and the paper says so.

## Baselines (M3)
1. Logistic regression (sklearn, default L2, standardised features) on: length, GRAVY (Kyte-Doolittle), net charge at pH 7 (Biopython),
   fraction charged (DEKR), aromatic fraction (FWY). Fit on the UESolDS TRAIN FASTA (never the test split).
2. SWI (published predictor of Gate 1a). 3. Majority class (AUC 0.5 by definition).
Winner statements only when 95% CIs (paired bootstrap of the difference) exclude 0; otherwise 'inconclusive'.

## Evaluation (M4)
Primary: AUROC + bootstrap CI. Secondary: AUPRC, MCC, accuracy, sensitivity, specificity at threshold 0.5 and at the Youden-optimal threshold
(optimal threshold chosen on the same data -> optimistic, flagged). Calibration: reliability curve (10 equal-frequency bins), Brier score,
ECE, calibration intercept/slope. The test set is exactly balanced (50% soluble); any prevalence mismatch with Aiki-Sol's training pool is a
confound for calibration and will be stated, with a prior-shift-corrected variant labelled exploratory.
Length: distribution + AUROC by length tercile.

## Probes (M5, all exploratory)
Tags: His6 and SUMO only. NOTE (found in source): Aiki-Sol documents stripping of His6/HiBit/Strep/linkers/initiator Met -- SUMO is NOT in
the documented list, so a SUMO shift is expected behaviour, not a bug; reported as such. His6 control: a composition-matched shuffle of
HHHHHH is HHHHHH itself, so the control is a length-matched random 6-mer without His runs (documented deviation). SUMO control: composition-
matched scrambled SUMO of identical length. Determinism (10 sequences twice), input robustness, nanobody subset (cut first if time-limited).

## Budget
API requests to the Aiki-Sol service capped at 1500 (every HTTP call logged in logs/api_requests.jsonl, including failures).
Wall-clock budget ~3 h from 2026-10-28 11:21 PT.
