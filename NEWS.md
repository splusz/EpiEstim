# EpiEstim 2.2-0

## MISC

* Added a `NEWS.md` file to track changes to the package. (#74, @zkamvar)
* Added tests for plotting with vdiffr. (#74, @zkamvar)
* Remove un-used dependencies that were added during hackout3: plyr, grid, 
  and plotly (#74, @zkamvar)
* Remove compare package from suggests and use test_that version. (#74, @zkamvar)
* Bump minimum required version of coarseDataTools to 0.6-4 (#71, @zkamvar)
* Incidence objects are now handled appropriately with accessors (#65, @zkamvar)

# EpiEstim 2.1-0

* TBD

# EpiEstim 2.0-0

* Changed function names to snake_case (only exception is that R remains capital letter to avoid confusion between the reproduction number R and the growth rate 
r) and to be more explicit; so `EstimateR` becomes `estimate_R`, `OverallInfectivity` becomes `oberall_infectivity`, `WT` becomes `wallinga_teunis`, and `DiscrSI` becomes `discr_si`. Names of arguments to these functions have also changed to snake_case. 
* Compatibility with `incidence` package: in the function `estimate_R`, the first argument, i.e. the incidence from which the reproduction number is calculated, can now be, either a vector of case counts (as in version 1.0-0) or an `incidence` object (see R package `incidence`).
* Accounting for imported cases: in the function `estimate_R`, the first argument, i.e. the incidence from which the reproduction number can now provide information about known imported cases: by specifying the first argument as either a dataframe with columns "local" and "imported", or an `incidence` object with two groups (local and imported, see R package `incidence`). This new feature is described in Thompson et al. Epidemics 2019 (currently in review).
* Additional methods available for function `estimate_R`: in addition to `non_parametric_si`, `parametric_si` and `uncertain_si`, which were already available in EpiEstim 1.0-0, two new methods have been added: `si_from_data` or `si_from_sample`. These allow feeding function `estimate_R` data on observed serial intervals (method `si_from_data`) or posterior samples of serial interval distributions obtained from such data (method `si_from_sample`). These new features are described in Thompson et al. Epidemics 2019 (currently in review).
* new argument `config` for `estimate_R` function: this is meant to minimise the number of arguments to function `estimate_R`; so arguments `method`, `t_start`, `t_end`, `n1`, `n2`, `mean_si`, `std_si`, `std_mean_si`, `min_mean_si`, `max_mean_si`, `std_std_si`, `min_std_si`, `max_std_si`, `si_distr`, `mean_prior`, `std_prior`, and `cv_posterior` are now specified as a group under this new `config` argument. Such a `config` argument must be of class `estimate_R_config` and can be obtained as a results of the new `make_config` function. 
 