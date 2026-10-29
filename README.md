# Independent external evaluation of Aiki-Sol (public beta)

A pre-registered, CPU-only check of Aikium's Aiki-Sol protein-solubility predictor. Written as an outside user's validation, not a review.

## Headline result
- **Primary metric (fixed before scoring):** AUROC of the marginal output on 500 sequences (250 soluble / 250 insoluble) sampled from the UESolDS independent test split = **0.858 (95% bootstrap CI 0.823-0.891)**.
- Baselines on the same sequences: logistic regression on 5 biophysical features 0.608 (0.557-0.654); published SWI 0.506 (0.453-0.555). Paired AUROC differences: +0.251 (0.194-0.308) and +0.353 (0.295-0.414).
- **Leakage caveat:** 324/500 (64.8%) evaluated sequences are exact MD5 matches to Aiki-Sol's released training list, so 0.858 is an **upper bound**. AUROC on the 176 non-matching sequences: 0.741 (0.668-0.807).
- Marginal output calibration: Brier 0.156, ECE 0.044, slope 0.85.
- Everything beyond the primary metric is exploratory. One benchmark, N = 500, mixed-protocol labels (the marginal output was chosen as primary for that reason).

## Layout
| Path | Contents |
|---|---|
| `PREREGISTRATION.md`, `.sha256` | Pre-specification written before any Aiki-Sol score was obtained (SHA-256 f804d4679f29b4bf2ec2aea2d2cad5bd6c7c06c10b319d5a758680e00251ea68) |
| `M0_feasibility.md` | Which access routes work (pip package does not exist on PyPI; REST and GitHub source do) and the output format |
| `paper/` | LaTeX source of the short write-up (`aikisol_eval_paper.tex`) |
| `figures/` | Fig. 1 ROC/PR with bootstrap bands, Fig. 2 reliability diagram |
| `slides/` | One slide with both figures |
| `code/` | All scripts. `lib.py` (loader, metrics, bootstrap, request-counting API wrapper) is the single shared module |
| `results/` | Raw predictions and result JSON/CSV |
| `data/` | Small derived inputs only: `controls.csv` (Gate 1b), `eval_sample.csv` (the 500 sequences), `zenodo_record.json` (file list used by the download scripts) |
| `logs/api_requests.jsonl` | Every HTTP call made to the Aiki-Sol service (126 of the 1,500-request cap) |

## Reproducing
Not included (fetch them; sources and licences are the owners'): 
1. UESolDS test/train FASTA from `github.com/Violet969/PLM_Sol` (`embedding_dataset/test_dataset.fasta` -> `data/plmsol_test_dataset.fasta`; `Train_dataset.fasta` -> `data/plmsol_train_dataset.fasta`).
2. Aiki-Sol source, commit `d8fcbf86` of `github.com/aikium-public/aiki-sol` -> `aikisol_src/` (the `aikisol/` package and `pyproject.toml`; unmodified).
3. Zenodo deposit 10.5281/zenodo.20151817 (record 20171586): `apache_tier__aikisol_v2_canonical_147k_full.pt` (2.6 GB, MD5 92ab0cc77dfaa78a66094dbf9493b38d) -> `data/zenodo/aikisol_v2_canonical_147k_full.pt`; `apache_tier__canonical_147k_md5.txt` -> `data/zenodo/`.  `code/download_checkpoint.py` and `code/fetch_more.py` do the Zenodo downloads.

Run order (from the repo root, so paths resolve): `gate1a.py` -> `gate1b_api.py` -> `gate1b_ci.py` -> `gate1b_supp_sample.py` + `local_aikisol.py` (Gate 1b supplementary) -> `make_eval_sample.py` -> `leakage_exact.py` -> `local_aikisol.py data/eval_sample.csv results/eval_local_predictions.csv` -> `evaluate.py` -> `make_figures.py` -> `probes_local.py` -> `probes_api.py` -> `probe_analysis.py` -> `build_slide.py`. Scoring 500 sequences took about 36 min on 8 CPU cores.

Python packages tested with: numpy 2.4, pandas, scipy 1.17, scikit-learn 1.9, biopython 1.88, torch 2.15, transformers 4.57, matplotlib, python-pptx (see `requirements.txt`). On Windows, installing torch into a deeply nested path fails on the 260-character limit; use a short install path.

## Deviations from the pre-registration (all disclosed in the paper)
- N = 500 sample of the 4,000-sequence split (CPU time).
- Gate 1a used SWI as the published predictor (SoluProt needs TMHMM/USEARCH; NetSolP is a 5.6 GB bundle).
- A supplementary direction check on 100 benchmark train-split sequences was added after the hand-picked Gate 1b controls were only weakly separated.
- The nanobody subset was not done (no public labelled VHH solubility set found).
- Compiled PDF: see `paper/` (the sandbox used could not run a TeX engine; compile `aikisol_eval_paper.tex` with any TeX Live / MiKTeX / Overleaf, with `figures/` next to `paper/`).
