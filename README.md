# Spectral CT for left atrial appendage thrombus: data and code

Data and code for the systematic review and Bayesian meta-analysis "Diagnostic performance of spectral CT for detection of left atrial appendage thrombus: a systematic review and meta-analysis". The protocol is registered on the Open
Science Framework: https://doi.org/10.17605/OSF.IO/3C2ZK. Methods and definitions are described in the article and in
its Online Resource 1.

## Contents

| Path | Content |
|---|---|
| `data/selection/records.csv` | The 90 unique records screened, with title and abstract decisions, full-text decisions, and the exclusion criterion (E1 to E6) |
| `data/extraction/extraction_items.csv` | Extracted items for the nine included reports: section, item, value, and location in the report |
| `data/analysis/dta_2x2.csv` | 2 × 2 tables entering each analysis, with their derivation and source in the report |
| `data/analysis/continuous.csv` | Spectral values by group used for the standardized mean differences |
| `data/analysis/quadas3.csv` | QUADAS-3 judgments and rationales for each estimate |
| `data/analysis/prisma.json` | Counts of the flow diagram |
| `search/search_strategies.md` | Full search strategies for the three databases |
| `code/12_dta_meta_analysis.R` | Bivariate models, sensitivity analyses, paired comparison, continuous values, consequences, and figures |
| `code/12b_posthoc_analyses.R` | Analyses labeled post hoc in the article |
| `results/` | Output tables of the two scripts |

Free-text notes in the extraction file are in English or Spanish. NR = not reported.

## Reproducing the analysis

R 4.5.2 with CmdStan 2.40.0 through cmdstanr 0.9.0, and the R packages writexl 1.5.4, jsonlite 2.0.0, systemfonts 1.3.2, ggrepel 0.9.8, ragg 1.5.2, tidybayes 3.0.7, distributional 0.9.0, ggdist 3.3.3, patchwork 1.3.2, ggplot2 4.0.1, tibble 3.3.0, purrr 1.2.0, readr 2.1.6, tidyr 1.3.2, dplyr 1.1.4, posterior 1.7.0, meta 8.5-0, mada 0.5.12, metafor 4.8-0, lme4 1.1-37, brms 2.23.0, plus
PRISMA2020, DiagrammeRsvg, rsvg, magick, and xml2 for the flow diagram. From the repository root:

```
Rscript code/12_dta_meta_analysis.R .
Rscript code/12b_posthoc_analyses.R .
```

Model fits are written to `fits/`, figures to `figures/`, and tables to `results/`. Both scripts fix the random seed.

## License

Code: MIT License (`LICENSE`). Data and documentation: Creative Commons Attribution 4.0 International (CC BY 4.0).

## Acknowledgment

This repository was prepared and uploaded with Claude (Anthropic).
