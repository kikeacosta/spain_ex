# From post-war hardship to longevity leader: a period and cohort analysis of the Spanish mortality paradox

Reproducible materials for the analysis of life expectancy convergence between
Spain and four reference countries (Sweden, Norway, the Netherlands and
Switzerland), 1908-2019, decomposing the convergence into period and cohort
contributions.

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
set of Figure S1 (Denmark, England and Wales, France).
Derived results that are small and not raw HMD data are kept in the repository:

| file | produced by | content |
|---|---|---|
| `data_inter/apc_cohort_effects.rds` | `05_fig3_cohort_effects_apc.R` | detrended cohort effects (rate ratios), five countries |
| `data_inter/tableS2_1_cohort_extrema.csv` | `05_fig3_cohort_effects_apc.R` | location and depth of the 1940s minimum and 1960s maximum |
| `data_inter/ap_residuals_ESP.rds` | `06_fig4_ap_residual_surfaces.R` | observed and model-fitted rates by age and year |
| `data_inter/cohort_conts_gap_changes.RData` | `07_fig5_cohort_contributions_gap.R` | cohort contributions to the change in the LE gap |
| `data_inter/simulation_scc.rds` | `08_figS3_simulation_scc.R` | simulated cohort contributions (three scenarios) |
| `data_inter/cohort_conts_gap_changes_altref.RData` | `09_figS1_alternative_reference.R` | same as Figure 5 with the alternative reference set |
| `data_inter/tableS4_convergence_1940_1960.csv` | `10_tableS4_convergence_1940_1960.R` | LE change and gap reduction 1940-1960 by sex and comparator |
| `data_inter/exceptional_cohorts.csv`, `disproportionate_cohorts.csv`, `exceptional_cohorts_simulation.csv` | `11_exceptional_cohorts.R` | contributions flagged by the exceptionality criterion, in the data and in the simulation (power check) |
| `data_inter/lexis_gap_contributions.rds` | `12_figS6_lexis_gap_contributions.R` | annual changes in the age-specific contributions to the e0 gap |

## Pipeline

Open `spain_ex.Rproj` (paths are relative to the project root) and run the
scripts in `code/` in order:

| script | manuscript output |
|---|---|
| `00_functions.R` | Arriaga decomposition and cohort-contribution functions (sourced by `07`) |
| `01_hmd_life_tables.R` | downloads e0 and life tables (inputs of Figures 1 and 5) |
| `02_hmd_deaths_exposures.R` | downloads deaths and exposures (inputs of Figures 3 and 4) |
| `03_fig1_life_expectancy_gap.R` | **Figure 1**: life expectancy 1908-2019 (A) and the gap between Spain and the reference mean (B) |
| `04_fig2_lexis_diagram.R` | **Figure 2**: Lexis diagram of the panel design for the cohort decomposition |
| `05_fig3_cohort_effects_apc.R` | **Figure 3** (Spain, ages 0-99), **Figure S2** (five countries, three age groups) and **Table S2.1** (depth of the 1940s minimum and 1960s maximum): detrended cohort effects, Age-drift-Period-Cohort models (`Epi::apc.fit`) |
| `06_fig4_ap_residual_surfaces.R` | **Figure 4**: Lexis surfaces of deviations from an Age-Period-drift model |
| `07_fig5_cohort_contributions_gap.R` | **Figure 5**: scaled cohort contributions to the change in the partial LE gap |
| `08_figS3_simulation_scc.R` | **Figure S3**: simulated cohort contributions under period-only and cohort-advantage scenarios, run through the same functions as Figure 5 |
| `09_figS1_alternative_reference.R` | **Figure S1**: Figure 5 with an alternative reference set (Denmark, England and Wales or the UK, France) |
| `10_tableS4_convergence_1940_1960.R` | **Table S4**: LE change and gap reduction in the 1940s, 1950s and 1940-1960, by sex and comparator |
| `11_exceptional_cohorts.R` | exceptionality criterion for cohort contributions (Methods, Step 3) and its power check on the simulation (Table S5.1); parameters at the top of the script |
| `12_figS6_lexis_gap_contributions.R` | **Figure S6**: Lexis surfaces of the annual changes in the age-specific contributions to the e0 gap, by sex and comparator |

Figures are written to `figures/manuscript/` and `figures/supplement/`.
Scripts `05`, `07` and `09` take several minutes (30 single-year APC factor
models; decompositions over 11 age groups x 3-4 countries x 80 years x 3
sexes, in both directions).

Numbers quoted in the text (life expectancy values, gaps and percentage
changes) can be checked with the objects printed by
`03_fig1_life_expectancy_gap.R`.

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

R (>= 4.3) with: `tidyverse`, `HMDHFDplus`, `Epi`, `MASS`, `splines`,
`MortalitySmooth` (`remotes::install_github("timriffe/MortalitySmooth")`),
`patchwork`. Exploratory scripts also use `cowplot`, `mgcv`, `scales`, `viridis`.

## Methods in brief

- **Figure 3.** Poisson APC models with single-year age, period and cohort as
  factors, "AdPC" parametrization (drift extracted with the unweighted inner
  product, `dr.extr = "1"` in `Epi`; period and cohort effects constrained to
  zero mean and zero slope), Spain 1940-2019, separately by sex for ages 0-49,
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
  1940-2019. The change in the contribution of age x from year y to y + 1 is
  assigned to the cohort aged x in year y and the changes are cumulated along
  each cohort's Lexis diagonal (only cohorts with a defined change at every
  age of the group), then divided by the largest absolute contribution within
  each sex, age group and country (scaled cohort contribution, SCC).
- **Figure S3.** Same functions applied to simulated mortality surfaces
  (Siler baseline declining 2% a year; Spanish excess log-mortality with a
  logistic catch-up centred in 1950, with or without a 15% advantage for the
  1945-1954 cohorts, or a linear convergence 1940-2000).
