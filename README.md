# From post-war hardship to longevity leader: a period and cohort analysis of the Spanish mortality paradox

Reproducible materials for the analysis of life expectancy convergence between
Spain and four reference countries (Sweden, Norway, the Netherlands and
Switzerland), 1940-2023 (life expectancy series shown from 1908), decomposing
the convergence into period and cohort contributions. The Supplementary
Material repeats the analysis with data ending in 2019, the period of the
original submission (S8).

## Data

All inputs come from the [Human Mortality Database](https://www.mortality.org)
(period life expectancy at birth, complete 1x1 life tables, and 1x1 death counts
and exposures). HMD data are **not** redistributed here. To reproduce the
analysis you need a free HMD account; store your credentials as environment
variables (e.g. in `~/.Renviron`):

```
hmd_us = "your_email"
hmd_pw = "your_password"
```

The download scripts write the raw data to `data_inter/` (ignored by git).
They fetch Spain, the four reference countries and the alternative reference
set of Supplementary Material S1 (Denmark, the United Kingdom, Canada and
Belgium; the United Kingdom series, HMD code `GBR_NP`, ends in 2022).
Derived results that are small and not raw HMD data are kept in the repository:

| file | produced by | content |
|---|---|---|
| `data_inter/apc_cohort_effects.rds` | `05_fig3_cohort_effects_apc.R` | detrended cohort effects (rate ratios), five countries |
| `data_inter/tableS2_1_cohort_extrema.csv` | `05_fig3_cohort_effects_apc.R` | location and depth of the 1940s minimum and 1960s maximum |
| `data_inter/ap_residuals_ESP.rds` | `06_fig4_ap_residual_surfaces.R` | observed and model-fitted rates by age and year |
| `data_inter/cohort_conts_gap_changes.RData` | `07_fig5_cohort_contributions_gap.R` | cohort contributions to the change in the LE gap |
| `data_inter/simulation_scc.rds` | `08_figS3_simulation_scc.R` | simulated cohort contributions (three scenarios) |
| `data_inter/cohort_conts_gap_changes_altref.RData`, `tableS1_1_altref_gap_1940_1960.csv` | `09_figS1_alternative_reference.R` | same as Figure 5 with the alternative reference set; gap with each alternative country in 1940, 1950 and 1960 (Table S1.1) |
| `data_inter/tableS4_convergence_1940_1960.csv` | `10_tableS4_convergence_1940_1960.R` | LE in Spain and the reference countries, gap and gap reduction 1940-1960 by sex and comparator |
| `data_inter/exceptional_cohorts.csv`, `disproportionate_cohorts.csv`, `block_pairs.csv`, `tableS5_1_exceptional_by_age_sex.csv`, `exceptional_cohorts_simulation.csv`, `exceptional_cohorts_sensitivity.csv`, `exceptional_cohorts_altref.csv` | `11_exceptional_cohorts.R` | exceptionality criterion: flagged contributions, replicated single peaks, opposite-signed pairs (block signature), counts by age interval, sex and sign (Table S5.1), power check on the simulation (Table S5.2), sensitivity (Table S5.3), flagged contributions with the alternative reference set |
| `data_inter/lexis_gap_contributions.rds` | `12_figS6_lexis_gap_contributions.R` | annual changes in the age-specific contributions to the e0 gap |
| `data_inter/tableS7_1_e0_vs_prewar_trend.csv` | `13_figS7_prewar_trend.R` | Spanish e0, its 1920-1935 linear trend and the deviation, selected years 1935-1960, by sex (Table S7.1) |
| `data_inter/apc_cohort_effects_2019.rds`, `tableS2_1_cohort_extrema_2019.csv`, `cohort_conts_gap_changes_2019.RData`, `cohort_conts_gap_changes_altref_2019.RData`, `tableS8_1_comparison_2019_2023.csv` | `14_figS8_prepandemic_2019.R` | the same objects estimated with data ending in 2019, and the comparison of the main indicators in both periods (Table S8.1) |

## Pipeline

Open `spain_ex.Rproj` (paths are relative to the project root) and run the
scripts in `code/` in order:

| script | manuscript output |
|---|---|
| `00_functions.R` | Arriaga decomposition and cohort-contribution functions (sourced by `07`, `08`, `09`, `12`, `14`), APC fitting and extrema functions (`05`, `14`), the exceptionality criterion (`11`, `14`) and the figure typeface (serif, `theme_bw()` and `theme_minimal()` redefined with `base_family = "serif"`; sourced by every figure script) |
| `01_hmd_life_tables.R` | downloads e0 and life tables (inputs of Figures 1, 5, S1.1, S7.1 and S8.2 and Tables S1.1, S4 and S7.1) |
| `02_hmd_deaths_exposures.R` | downloads deaths and exposures (inputs of Figures 3, 4, S2 and S8.1) |
| `03_fig1_life_expectancy_gap.R` | **Figure 1**: life expectancy 1908-2023 (A) and the gap between Spain and the reference mean (B) |
| `04_fig2_lexis_diagram.R` | **Figure 2**: Lexis diagram of the panel design for the cohort decomposition |
| `05_fig3_cohort_effects_apc.R` | **Figure 3** (Spain, ages 0-99), **Figures S2.1-S2.5** (reference countries and the five countries overlaid, ages 0-99, 0-49 and 50-99) and **Table S2.1** (depth of the 1940s minimum and 1960s maximum): detrended cohort effects, Age-drift-Period-Cohort models (`Epi::apc.fit`) |
| `06_fig4_ap_residual_surfaces.R` | **Figure 4**: Lexis surfaces of deviations from an Age-Period-drift model |
| `07_fig5_cohort_contributions_gap.R` | **Figure 5**: scaled cohort contributions to the change in the partial LE gap |
| `08_figS3_simulation_scc.R` | **Figure S3.1**: simulated cohort contributions under period-only and cohort-advantage scenarios (1940-2019), run through the same functions as Figure 5 |
| `09_figS1_alternative_reference.R` | **Figure S1.1** and **Table S1.1**: Figure 5 with an alternative reference set (Denmark, the United Kingdom to 2022, Canada and Belgium), and the gap with each of them in 1940, 1950 and 1960 |
| `10_tableS4_convergence_1940_1960.R` | **Table S4**: LE in Spain and the reference countries, gap and gap reduction in the 1940s, 1950s and 1940-1960, by sex and comparator |
| `11_exceptional_cohorts.R` | exceptionality criterion for cohort contributions (Methods, Step 3; **Tables S5.1-S5.3**): robust score against the ten neighbouring cohorts, two-sided, edges excluded, replication across comparators and age intervals, block-pair signature, power check on the simulation and sensitivity; also applied to the alternative reference set (S1); parameters at the top of the script |
| `12_figS6_lexis_gap_contributions.R` | **Figure S6.1**: Lexis surfaces of the annual changes in the age-specific contributions to the e0 gap, by sex and comparator |
| `13_figS7_prewar_trend.R` | **Figure S7.1** and **Table S7.1**: Spanish e0 against a linear trend fitted to 1920-1935 and extrapolated to 1960, by sex |
| `14_figS8_prepandemic_2019.R` | **Figures S8.1-S8.2** and **Table S8.1**: results with data ending in 2019 (detrended cohort curves for Spain in both periods, Figure 5 to 2019, and the main indicators in both periods) |

Figures are written to `figures/manuscript/` and `figures/supplement/`.
Scripts `05`, `07`, `09` and `14` take several minutes to half an hour each
(30 single-year APC factor models in `05` and again in `14`; decompositions
over 11 age groups x 4 countries x 84 years x 3 sexes, in both directions, in
`07`, `09` and twice in `14`).

Numbers quoted in the text (life expectancy values, gaps and percentage
changes) can be checked with the objects printed by
`03_fig1_life_expectancy_gap.R` and `10_tableS4_convergence_1940_1960.R`.

### Exploratory material

`code/exploratory/` holds scripts that are not needed for the manuscript:
descriptive plots, Lexis surfaces of mortality change, earlier versions of the
decomposition, and additional views of the model deviations. They source
`code/exploratory/00_functions_exploratory.R` (which in turn sources the
manuscript functions). `code/exploratory/previous_versions/` keeps the scripts
as they were before this clean-up, and `code/exploratory/deprecated/` older
developments. Their outputs go to `figures/exploratory/` (ignored by git).
`code/exploratory/python/simulation_SCC.py` is the original Python version of
the simulation (numpy + matplotlib), kept for reference; `08_figS3_simulation_scc.R`
is its R port.

## Requirements

R (>= 4.3) with: `tidyverse` (`ggplot2` >= 4.0, so that text geoms inherit the
theme typeface), `HMDHFDplus`, `Epi`, `MASS`, `splines`,
`MortalitySmooth` (`remotes::install_github("timriffe/MortalitySmooth")`),
`patchwork`. Exploratory scripts also use `cowplot`, `mgcv`, `scales`, `viridis`.

## Methods in brief

- **Figure 3.** Poisson APC models with single-year age, period and cohort as
  factors, "AdPC" parametrization (drift extracted with the unweighted inner
  product, `dr.extr = "1"` in `Epi`; period and cohort effects constrained to
  zero mean and zero slope), 1940-2023, separately by sex for ages 0-49,
  50-99 and 0-99, for Spain and the four reference countries. The choice of
  inner product only tilts the detrended curves linearly (second differences
  are invariant), so cohort effects are read through their curvature; Table
  S2.1 measures the depth of each extremum relative to the line joining the
  neighbouring turning points (prominence in log RR), which is invariant to
  that tilt.
- **Figure 4.** `curvature_dAPC()` (an adaptation of `Epi::apc.fit`) fits the
  full APC model and refits it with the detrended cohort effects set to zero,
  leaving an Age-drift-Period model; the surfaces show the relative deviation
  of observed rates from that model.
- **Figure 5.** For each 10-year age group, the gap in partial life expectancy
  (reference country minus Spain) is decomposed by single age with the
  symmetric Arriaga method (average of the two directions) every year
  1940-2023. The change in the contribution of age x from year y to y + 1 is
  assigned to the cohort aged x in year y and the changes are cumulated along
  each cohort's Lexis diagonal (only cohorts with a defined change at every
  age of the group, i.e. cohorts 1940 - a to 2013 - a for the group starting
  at age a), then divided by the largest absolute contribution within each
  sex, age group and country (scaled cohort contribution, SCC). The
  alternative reference set (Figure S1.1) uses the same window, except that
  the comparison with the United Kingdom ends in 2022.
- **Exceptionality criterion (script 11, Tables S5.1-S5.3).** A cohort
  contribution is exceptional when its robust z-score against the ten
  neighbouring cohorts (median and scaled MAD, same sex, age interval and
  comparator) exceeds 3 in absolute value; the first and last five cohorts of
  each interval are not tested. A single peak is disproportionate when the
  same sign recurs in at least three comparators within an interval and in at
  least two intervals (period changes fall on different cohorts in different
  intervals, cohort deviations on the same one). A block of advantaged
  cohorts shows as a narrowing peak at its entry and a widening peak at its
  exit, so opposite-signed pairs 2-15 cohorts apart are searched for and
  replicated in the same way. The criterion detects discontinuities, not
  gradual cohort advantages. It is validated on the simulation (scenario B:
  edges 1944 and 1954) and reported for seven settings of scale, window and
  threshold. The criterion and its parameters are the same for 1940-2023 and
  for 1940-2019 (S8).
- **Figure S3.** Same functions applied to simulated mortality surfaces
  (Siler baseline declining 2% a year; Spanish excess log-mortality with a
  logistic catch-up centred in 1950, with or without a 15% advantage for the
  1945-1954 cohorts, or a linear convergence 1940-2000), decomposed over
  1940-2019 as in the original submission.
- **Figure S7.** A linear trend of Spanish life expectancy at birth is fitted
  by ordinary least squares to 1920-1935, separately for both sexes, women and
  men, and extrapolated to 1960; the deviation is observed minus trend.
- **S8 (script 14).** Every analysis is re-estimated with data ending in 2019
  (APC models, decompositions with both reference sets, exceptionality
  criterion) with the same functions and parameters; the 1940-2023 results
  are read from the outputs of scripts 05, 07 and 09.
