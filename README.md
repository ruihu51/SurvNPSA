# SurvNPSA

`SurvNPSA` implements a nonparametric causal sensitivity analysis framework for time-to-event data. It provides  
- nonparametric bounds and valid inference for survival contrasts, including the difference in survival curves and the difference in restricted mean survival times (RMST);
- summary metrics to help practitioners understand the extent of unobserved confounding required to explain away the observed causal effect

## Installation

You can install the development version of **SurvNPSA** from GitHub:

```r
# install.packages("pak")
pak::pak("ruihu51/SurvNPSA")
```

## Choosing Time Points

Most users can start without manually choosing time points:

```r
np_out <- np_surv(
  time = time,
  event = event,
  treat = treat,
  confounders = confounders
)
```

When `fit.times` is not supplied, `SurvNPSA` starts from the CFsurvival-style
idea of using positive observed follow-up times before the largest observed
event time. For large datasets, it then keeps a compact set of representative
times and avoids the very late tail when overall censoring support is weak. The
censoring-support cutoff is a practical safety default to reduce unstable tail
estimates, not a change to the statistical formulas.

The main time options are:

- `fit.times`: the main analysis grid where survival and effect estimates are computed.
- `nuisance.options$eval.times`: the internal prediction grid for nuisance survival models; most users should not set this.
- `bound.options$plot.times`: times reported for sensitivity bounds; defaults to all fitted times.
- `rv.options$rv.times`: times where RV/MIRV are computed; defaults to about five representative fitted times.
- `uniform.cutpoint`: defines the uniform inference window; it is not a set of pointwise output times.
- `rmst.options$fit.times.rmst`: RMST horizons; choose these only when running `npsa_surv(rmst = TRUE)`.

## Codes Structure

```
R/
├── npsa_main.R             
├── estimation/             
│   ├── estimate_nuisances.R
│   ├── estimate_obs_components.R
│   ├── estimate_rmst.R
├── senspar/                
│   ├── npsa_senspar.R
├── sensitivity_utils/          
│   ├── npsa_summary.R
│   ├── npsa_utils.R
```
