# SampleMultiCat

Bilingual (English/Spanish) R package for sample size calculation in regression models with polytomous outcomes.

## Features

- **`SampleOrdinal()`**: Sample size for proportional odds ordinal logistic regression (worst cumulative cut EPP approach). Optional power-based calculation via `Hmisc::posamsize`.
- **`SampleMultinomial()`**: Sample size for multinomial logistic regression (EPPm criterion from de Jong et al. 2019).
- **`ValidateSample()`**: Post-hoc validation of EPP/EPPm for published studies.
- Scenario tables (EPP/EPPm = 10, 20, 30, 40, 50) for planning.
- Bilingual output: `language = 'es'` (Spanish) or `language = 'en'` (English).

## Important: what is `k`?

The parameter `k` refers to **total regression parameters** (degrees of freedom consumed by all predictors), **not** the number of variables.

| Variable type | How to count |
|---|---|
| Continuous | 1 parameter per variable |
| Categorical (m levels) | m - 1 parameters |

**Example:** A model with 3 continuous variables + 1 categorical variable with 4 levels:
`k = 3 + (4-1) = 6`

This follows the EPP definitions from Peduzzi et al. (1996), Harrell (2015), and de Jong et al. (2019), who base the criterion on regression parameters (degrees of freedom).

## Installation

```r
# install.packages('devtools')
devtools::install_github('VicVePo/SampleMultiCat')
```

## Quick examples

```r
library(SampleMultiCat)

# Ordinal: depression scale (none, mild, moderate, severe)
# Model: 5 continuous + 1 categorical (3 levels) = 5 + 2 = 7 parameters
SampleOrdinal(k = 7, probs = c(0.50, 0.20, 0.20, 0.10), EPP = 20)

# Ordinal with power component
SampleOrdinal(k = 7, probs = c(0.50, 0.20, 0.20, 0.10), EPP = 20, OR_min = 1.5)

# Multinomial: 5 tumor subtypes
# Model: 10 continuous + 1 categorical (4 levels) = 10 + 3 = 13 parameters
SampleMultinomial(k = 13, probs = c(0.40, 0.25, 0.20, 0.10, 0.05), EPPm = 20)

# Scenario table
SampleOrdinal(k = 7, probs = c(0.50, 0.20, 0.20, 0.10), scenarios = TRUE)

# Post-hoc validation
ValidateSample(observed_counts = c(500, 200, 150, 50), k = 8, model = 'ordinal')
```

## Changelog

### v0.5.0
- Fixed floating-point precision bug (v0.4.0): `ceiling()` on imprecise divisions (e.g. `60 / 0.0999...`) produced N+1 instead of N. Now uses `safe_ceiling()` = `ceiling(round(x, 8))`.
- Changed k label from Predictors to Parameters throughout (v0.4.1).
- Unified category symbol to (J) in both models; renamed EPV to EPP and EPVm to EPPm for consistency with k = parameters (v0.5.0).

## References

- Peduzzi P, et al. (1996). J Clin Epidemiol 49:1373-1379.
- Harrell FE (2015). Regression Modeling Strategies. Springer.
- de Jong VMT, et al. (2019). Stat Med 38:1601-1619.

## Author

Victor J. Vera-Ponce ([ORCID](https://orcid.org/0000-0003-4075-9049))
