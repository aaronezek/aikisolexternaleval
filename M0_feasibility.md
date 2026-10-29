# Milestone 0 - feasibility of access routes (written 2026-10-28T12:49:19-07:00)

| Route | Result | Evidence |
|---|---|---|
| `pip install aikisol` | **Does not work as advertised**: package is not on PyPI (HTTP 404 on pypi.org/pypi/aikisol/json). | README and landing page both say `pip install aikisol`. |
| Install from GitHub source | **Works** (no pip needed): 4 source files of github.com/aikium-public/aiki-sol (commit d8fcbf86) run unmodified on CPU/Windows with torch 2.15 + transformers 4.57. | `code/local_aikisol.py`; 40 sequences scored in 67 s after a 14 s model load. |
| Model weights | Zenodo 10.5281/zenodo.20151817 (record 20171586), `aikisol_v2_canonical_147k_full.pt`, 2,605,701,359 bytes, MD5 92ab0cc77dfaa78a66094dbf9493b38d verified against Zenodo checksum. Needs zenodo.org + huggingface.co (tokenizer/config) access. | `logs/download_ckpt.log` |
| Docker (`ghcr.io/aikium-public/aiki-sol:full`) | **Not tried**: Docker is not installed on this machine and ghcr.io is not on the sandbox allowlist. | `where docker` returns nothing |
| REST `GET /lookup_protein?seq=` | **Works** | 44/44 HTTP 200, median latency 1.3 s (p95 2.8 s, max 3.0 s) at 1 s spacing; no throttling observed. No burst test was run (deliberately not hammering the server). |
| REST `POST /predict_fasta` (json / csv) | **Works** (tested with 2 sequences each). | `logs/api_requests.jsonl` |
| Web UI | Not used (and not scraped). | - |

**API vs local package:** on the 40 Gate-1b controls, all 7 score columns agree to <= 1.85e-6 absolute and `normalized_sequence` is identical, so the hosted service and the Zenodo canonical checkpoint are the same model for practical purposes. Bulk scoring therefore uses the local route (no API budget used); the API is used for the gates and probes that test the hosted tool itself.

## Output format (documented + observed)
`/lookup_protein` returns JSON: `mean_prob` (mean of 5 regime outputs), `prob_marginal` (dedicated 6th head), `per_stringency` = {`3000g_10min`, `6000g`, `32000g`, `eSol_21600g_30min`, `100000g`}, `normalized_sequence` (what was actually scored; tags stripped, initiator Met re-added).
`/predict_fasta` returns `{"n":..,"predictions":[{"id","mean_prob","prob_marginal","normalized_sequence","prob_<regime>"...}]}` (or the same columns as CSV).
Range: all outputs are sigmoid probabilities in [0, 1]. Meaning per README/landing page: `prob_<regime>` = P(soluble | that centrifugation protocol); `prob_marginal` = P(soluble | sequence) when the protocol is unspecified ("the recommended single number"); `mean_prob` = unweighted mean of the five regime outputs. Higher = more soluble. Landing page and code both describe outputs as "calibrated".
Documented limits: sequences truncated to 1,022 residues; tags (His6, HiBit, Strep, common linkers, initiator Met) stripped by `normalize_sequence` (SUMO is NOT in the documented list); batch preflight raises when >5% of a batch retain His6/HiBit; predictor is for E. coli; not for yield/aggregation kinetics/stability/activity.

## API request ledger
44 requests so far (4 format-documentation calls in `m0_explore.py`, 40 Gate-1b calls). Cap = 1500. Source of truth: `logs/api_requests.jsonl` (one line per HTTP call, including failures).
