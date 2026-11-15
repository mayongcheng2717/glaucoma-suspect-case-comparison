# Glaucoma suspects carry about half the genetic liability of confirmed disease

Code, derived tables and figures for:

> Wang X, Li Y, Ma Y, Yang C, Zheng Y, Gan H, Yang X. *Glaucoma suspects carry about half the genetic
> liability of confirmed disease: a genome-wide case–case comparison.* (manuscript under review)

A genome-wide **case–case** comparison of 16,604 glaucoma suspects against 11,637 confirmed
primary open-angle glaucoma cases in FinnGen release 13. The contrast uses the alt-allele
frequency **among cases**, which FinnGen reports separately for every endpoint, so it requires
no control data. Two independent loci separate the two groups (9p21/CDKN2B, 1q24/TMCO1); both
behave as a continuum rather than as stage markers, and both sit at unity in wet and dry
age-related macular degeneration, which serve as internal negative controls.

## Layout

```
code/            analysis scripts (see Pipeline below)
code/dse/        shared package: FinnGen access, remote tabix reader, statistics, plotting
code/novelty/    prior-art / novelty checks run before submission
tables/          derived tables reported in the manuscript
figures/         Figures 1-7 (PNG and PDF)
manuscript/      manuscript PDF
```

**Repository:** <https://github.com/mayongcheng2717/glaucoma-suspect-case-comparison>

## Data sources

| Data | Access |
| --- | --- |
| FinnGen release 13 summary statistics | https://r13.finngen.fi ; public bucket `https://storage.googleapis.com/finngen-public-data-r13/summary_stats/` |
| GTEx v8 cis-eQTL summary statistics | eQTL Catalogue, https://www.ebi.ac.uk/eqtl |
| Ensembl / GWAS Catalog / Europe PMC | public REST APIs, queried by `code/dse/` |

No individual-level data were accessed. Only publicly released summary statistics are used.

## Pipeline

Run in this order (each step writes into `results/gwcas/`; `extract_casearm.py` fills
`cache/casearm/` first):

| Step | Script | Produces |
| --- | --- | --- |
| 1 | `code/extract_casearm.py` | per-endpoint case-arm allele frequencies (streamed, nothing stored beyond what is needed) |
| 2 | `code/gwcas.py` | genome-wide case–case scan, lead variants, λ, top-200k scan output |
| 3 | `code/qq_proper.py` | QQ data for the calibration figure |
| 4 | `code/lead_specificity.py` | the two lead variants across every ocular endpoint |
| 5 | `code/rebuild_or_table.py` | case-control odds-ratio table (Table 1) |
| 6 | `code/eqtl_anchor.py` | cis-eQTL anchoring at both loci |
| 7 | `code/glauc_figs_final.py` | Figures 1–7 |

## Local paths

The scripts were written inside a working tree rooted at `D:\eye100` and reference it directly
(`CACHE`, `OUT`, `RES`, `DEL` constants, and `ROOT` in `code/dse/finngen.py`). **Change that one
root to your own path before running**, and download the endpoints you need into
`<root>/cache/casearm/` in step 1. Nothing else is machine-specific.

## Not included

The four top-200k scan files (`.tsv.gz`, ~34 MB total) that back Figure 1 and the calibration
contrasts ship with the submission as supplementary files rather than in this repository.
They can be regenerated with `code/gwcas.py`.

## Requirements

Python 3.12 with numpy, pandas, scipy, matplotlib and requests (see `requirements.txt`).
No R installation is needed: every statistic (logit variance, BH-FDR, genomic control) is
implemented from first principles in `code/dse/stats.py`.

## License

CC BY 4.0 — see `LICENSE`.

## Citation

If you use this code, please cite the manuscript above. A DOI for this repository will be
minted on the first release (Zenodo).
