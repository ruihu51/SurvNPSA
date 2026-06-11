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
time.out <- npsa_times(time, event)
summary(time.out)
plot(time.out)

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
- `rv.options$rv.times`: times where RV/MIRV are computed; defaults to `plot.times` when users set `plot.times`, otherwise about five representative fitted times.
- `bound.options$uniform.cutpoint` or `bound.options$uniform.window`: defines the uniform inference window; it is not a set of pointwise output times.
- `rv.options$uniform.window` or `rv.options$uniform.cutpoint`: optional URV-only window override.
- `rmst.options$fit.times.rmst`: RMST horizons; choose these only when running `npsa_surv(rmst = TRUE)`.

Common time-setting scenarios:

```r
# 1. Use all defaults
time.out <- npsa_times(time, event)

np_out <- np_surv(time, event, treat, confounders)
npsa_out <- npsa_surv(time, event, treat, confounders,
                      result = np_out$result)

# Expected:
# fit.times  = automatic, about 50 representative times
# eval.times = automatic, about 200 internal times
# plot.times = all fit.times
# rv.times   = about 5 representative fit.times
```

```r
# 2. User specifies analysis, plot, and RV times
my_fit_times <- seq(0.1, 1.2, by = 0.1)
my_plot_times <- c(0.2, 0.6, 1.0)
my_rv_times <- c(0.2, 0.6, 1.0)

time.out <- npsa_times(time, event,
                       fit.times = my_fit_times,
                       plot.times = my_plot_times,
                       rv.times = my_rv_times)

np_out <- np_surv(time, event, treat, confounders,
                  fit.times = my_fit_times,
                  np.options = list(plot.times = my_plot_times))
npsa_out <- npsa_surv(time, event, treat, confounders,
                      fit.times = my_fit_times,
                      result = np_out$result,
                      bound.options = list(plot.times = my_plot_times),
                      rv.options = list(rv.times = my_rv_times))

# Expected:
# fit.times  = seq(0.1, 1.2, by = 0.1)
# eval.times = automatic, covers 0 to 1.2
# plot.times = c(0.2, 0.6, 1.0)
# rv.times   = c(0.2, 0.6, 1.0)
```

```r
# 3. User specifies only the main analysis times
my_fit_times <- c(0.2, 0.4, 0.6, 0.8)
time.out <- npsa_times(time, event, fit.times = my_fit_times)

np_out <- np_surv(time, event, treat, confounders,
                  fit.times = my_fit_times)
npsa_out <- npsa_surv(time, event, treat, confounders,
                      fit.times = my_fit_times,
                      result = np_out$result)

# Expected:
# fit.times  = c(0.2, 0.4, 0.6, 0.8)
# eval.times = automatic, covers 0 to 0.8
# plot.times = c(0.2, 0.4, 0.6, 0.8)
# rv.times   = c(0.2, 0.4, 0.6, 0.8)
```

In general, choose `fit.times` first. Then choose `plot.times` and
`rv.times` from `fit.times`.

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
