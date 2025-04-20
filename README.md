# SurvNPSA

`SurvNPSA` implements a nonparametric causal sensitivity analysis framework for time-to-event data. Specifically, it provides   
- nonparametric bounds and valid inference for survival contrasts, including the difference in survival curves and the difference in restricted mean survival times (RMST);
- summary metrics to help practitioners understand the extent of unobserved confounding required to explain away the observed causal effect

## Installation

You can install the development version of **SurvNPSA** from GitHub:

```r
# install.packages("devtools")
devtools::install_github("ruihu51/SurvNPSA")
```

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

