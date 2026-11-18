# Compound drought–wildfire indicators: analysis code

R code for the manuscript "An Empirical Framework for Identifying Indicators of Compound Drought–Wildfire Risk". The scripts link Korea Forest Service wildfire damage records (2000–2025) to antecedent drought indices (SPI, SPEI) and fire weather from KMA ASOS stations, and estimate how drought conditions relate to burned area, suppression duration and suppression resources at the incident level.

## Structure

| File | Purpose |
|---|---|
| `R/run_all.R` | Runs the full pipeline from the repository root |
| `R/00_download.R` | KMA ASOS daily data and station metadata via the KMA API Hub |
| `R/download_asos.R` | Example session: ASOS download, station remapping, panel build and checks |
| `R/01_spi.R` | SPI (gamma with zero-precipitation mixture) and SPEI (Hargreaves PET, log-logistic) |
| `R/10_sigungu.R` | Municipality (si/gun/gu) name harmonisation |
| `R/11_parse_damage_ledger.R` | Parses the wildfire damage ledger; merges duplicate entries of the same fire complex |
| `R/12_station_remap.R` | Station audit and municipality–station matching by record length |
| `R/13_build_panel.R` | Incident-level panel: drought indices for the month before ignition, fire weather, suppression resources |
| `R/14_analysis_q2.R` | Fixed-effects regressions, quantile regression, propensity score matching, resource and duration models, sensitivity grid |
| `R/16_make_tables_figures.R` | Tables and quantile-regression results |
| `R/17_figures_revision_rev.r` | Estimates behind Figures 2–5 and S1 |
| `R/figs/Figure3_redraw.R`, `Figure5_redraw.R`, `FigureS1_redraw.R` | Final versions of Figures 3, 5 and S1 |
| `sigungu_centroid.csv` | Municipality centroids used for station matching |

## Data

Raw data are not included. Place the following in `data_raw/`:
Author can provide input data if you request
- Wildfire damage ledger 2000–2025 (`산불피해대장`, Korea Forest Service)
- Suppression resource history (`진화자원이력.xlsx`, Korea Forest Service)
- KMA ASOS daily observations and station metadata (`kma_asos_daily.csv`, `station_meta.csv`), downloadable with `R/00_download.R`

The KMA API Hub key is read from the environment. Copy `.Renviron.example` to `.Renviron` and add your key.

## Running

From the repository root:

```r
source("R/run_all.R")
```

Outputs are written to `data_out/` (panel, tables, report) and `figs/revision_final/` (figures).

Required packages: data.table, sandwich, lmtest, quantreg, MASS, readxl, tidyverse, patchwork, scales.

## License

MIT License. See `LICENSE`.
