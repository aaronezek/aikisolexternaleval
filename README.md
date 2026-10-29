# Independent external evaluation of Aiki-Sol (public beta)

A pre-registered, CPU-only check of Aikium's Aiki-Sol protein-solubility predictor, written as an outside user's validation rather than a review. The short paper source is in `paper/`; a one-slide summary is in `slides/`.

## Headline result
- **Primary metric (fixed before scoring):** AUROC of Aiki-Sol's marginal output on 500 sequences (250 soluble, 250 insoluble) sampled from the UESolDS independent test split: **0.858 (95% bootstrap CI 0.823-0.891)**.
- Same sequences, baselines: logistic regression on 5 biophysical features 0.608 (0.557-0.654); published SWI 0.506 (0.453-0.555). Paired AUROC differences: +0.251 (0.194-0.308) and +0.353 (0.295-0.414).
- **Leakage caveat:** 324 of 500 evaluated sequences (64.8%) are exact MD5 matches to Aiki-Sol's released training list, so 0.858 is an upper bound. AUROC on the 176 non-matching sequences is 0.741 (0.668-0.807).
- The marginal output was reasonably calibrated: Brier 0.156, ECE 0.044, calibration slope 0.85.
- Everything beyond the primary metric is exploratory. One benchmark, N = 500, mixed-protocol labels (which is why the marginal output was chosen as primary).

## Repository layout
| Path | Contents |
|---|---|
| `PREREGISTRATION.md`, `PREREGISTRATION.sha256` | Analysis plan written before any Aiki-Sol score was obtained; `sha256sum PREREGISTRATION.md` reproduces the recorded hash (f804d467...) |
| `M0_feasibility.md` | Which access routes worked and the observed output format |
| `paper/` | LaTeX source (`aikisol_eval_paper.tex`); build with any TeX distribution or Overleaf, with `figures/` available as `../figures/` |
| `figures/` | Fig. 1 (ROC and precision-recall with bootstrap bands), Fig. 2 (reliability diagram) |
| `slides/` | One-slide summary of both figures |
| `code/` | All scripts; `lib.py` holds the shared loaders, metrics, bootstrap and the request-counting API wrapper |
| `results/` | Raw predictions and result tables |
| `data/` | `controls.csv`, `eval_sample.csv` (the 500 sequences), `zenodo_record.json` (file list used by the download scripts) |
| `logs/api_requests.jsonl` | Every HTTP call made to the hosted service (126 of the 1,500-request cap) |

## Reproducing
Fetch the inputs that are not included:
1. UESolDS FASTA files from `github.com/Violet969/PLM_Sol` (`embedding_dataset/test_dataset.fasta` to `data/plmsol_test_dataset.fasta`, and `Train_dataset.fasta` to `data/plmsol_train_dataset.fasta`).
2. Aiki-Sol source at commit `d8fcbf86` of `github.com/aikium-public/aiki-sol` into `aikisol_src/` (the `aikisol/` package, unmodified).
3. The Zenodo deposit 10.5281/zenodo.20151817 (record 20171586): the canonical checkpoint (2.6 GB, MD5 92ab0cc77dfaa78a66094dbf9493b38d) saved as `data/zenodo/aikisol_v2_canonical_147k_full.pt`, and `apache_tier__canonical_147k_md5.txt` in `data/zenodo/`. `code/download_checkpoint.py` and `code/fetch_more.py` fetch these.

Run from the repository root: `gate1a.py`, `gate1b_api.py`, `gate1b_ci.py`, `gate1b_supp_sample.py` with `local_aikisol.py`, `make_eval_sample.py`, `leakage_exact.py`, then `local_aikisol.py data/eval_sample.csv results/eval_local_predictions.csv`, `evaluate.py`, `make_figures.py`, `probes_local.py`, `probes_api.py`, `probe_analysis.py`, `build_slide.py`. Scoring 500 sequences took about 36 minutes on 8 CPU cores. `requirements.txt` lists the packages used.

## Deviations from the analysis plan
- N = 500 sampled from the 4,000-sequence split, because scoring was CPU-only.
- The published predictor re-run for the pipeline check was SWI; SoluProt needs TMHMM and USEARCH, and NetSolP is a 5.6 GB bundle.
- A supplementary direction check on 100 benchmark train-split sequences was added after the 40 hand-picked controls were only weakly separated.
- No nanobody subset: no public labelled VHH solubility set was found.

## Attribution and licence
Code in this repository is released under the MIT licence (see `LICENSE`). The sequences in `data/eval_sample.csv` and the derived result files come from the public UESolDS benchmark (Zhang et al., Briefings in Bioinformatics 25(5):bbae404, 2024; repository `Violet969/PLM_Sol`). Aiki-Sol code and weights belong to Aikium and are not included here.
