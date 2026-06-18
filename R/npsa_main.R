#' Nonparametric Sensitivity Analysis for Survival Outcomes
#'
#' This function estimates the treatment effect on survival outcomes with
#' nonparametric sensitivity analysis under unmeasured confounding, optionally
#' including restricted mean survival time (RMST) analysis and robustness value computation.
#'
#' If \code{fit.times} is not supplied, the function chooses a compact analysis
#' grid from positive observed follow-up times before the largest observed event
#' time, following the CFsurvival default idea. For large data, this grid is
#' reduced to representative times and capped before the tail where overall
#' censoring support is very weak. The censoring-support cutoff is a practical
#' safety default, not a change to the statistical formulas. Most users do not
#' need to set \code{nuisance.options$eval.times}; it is an internal prediction
#' grid and is chosen automatically when omitted.
#'
#' Sensitivity bound curves are computed over \code{fit.times} for plotting.
#' \code{bound.options$report.times} controls reported time summaries and the
#' default \code{rv.options$rv.times} when RV/MIRV times are not supplied.
#'
#' @param time Numeric vector of event or censoring times.
#' @param event Numeric vector of event indicators (1 = event, 0 = censored).
#' @param treat Numeric vector of treatment assignment indicators (1 = treated, 0 = control).
#' @param confounders Matrix or data frame of observed confounders (observed covariates).
#' @param fit.times Optional numeric vector of times at which the survival
#'   contrasts are estimated. If \code{NULL}, a compact default grid is chosen
#'   from the observed follow-up times.
#' @param nuisance.options List of options for nuisance estimation.
#' @param target.options List of options for target parameter estimation.
#'   May include \code{psi.type} and \code{tau.type}.
#' @param bound.options List of options for reporting pointwise and uniform bounds.
#'   The \code{transform} option is also used for pointwise MIRV calculation.
#'   May include \code{report.times}, \code{uniform.cutpoint}, or an exact
#'   \code{uniform.window}. Bounds are still computed over \code{fit.times}
#'   for plotting.
#' @param rv.options List of options for robustness value computation. May include
#'   \code{rv.times}, \code{rho}, and \code{theta}. The \code{rho} option is
#'   also used for sensitivity-bound calculations so bounds and RV use the same
#'   correlation setting. The \code{uniform.window} and \code{uniform.cutpoint}
#'   options can be used here for a URV-specific window; otherwise URV uses the
#'   same window as \code{bound.options}.
#'   If \code{rv.times = NULL}, RV/MIRV are computed at about five
#'   representative fitted times.
#' @param rmst Logical; if TRUE, estimate RMST and its bounds inference as well.
#' @param rmst.options List of options for RMST estimation. If \code{rmst = TRUE}
#'   and \code{fit.times.rmst} is not supplied, the default RMST horizon is
#'   \code{max(fit.times)} with a message reminding users that RMST horizons
#'   need analyst interpretation.
#' @param sens.options List of options for sensitivity parameter simulation.
#'   Use either \code{pct_drop} for percentage-based dropping or \code{num_drop}
#'   for exact numbers of covariates to drop, but not both. The default uses
#'   \code{pct_drop = c(0.3, 0.7)}. The simulation also includes the benchmark
#'   drop sizes \code{1} and \code{ceiling(0.5 * n_var)} when possible.
#'   Use \code{senspar.save.path} to save generated sensitivity parameters to
#'   a custom path. Use \code{senspar.only = TRUE} to stop after sensitivity
#'   parameter simulation.
#' @param result Optional precomputed result object (e.g., containing nuisances).
#' @param var_names Character vector of confounder variable names.
#' @param verbose Logical; if TRUE, print system timestamps for each estimation step.
#' @param save Logical; if TRUE, save intermediate results.
#'
#' @return A list of class \code{npsa_surv} containing:
#' \describe{
#'   \item{result}{Estimated observable components and IFs, such as observed survival differences and rmst differences.}
#'   \item{bounds.df}{Estimated bounds on survival contrasts over time.}
#'   \item{senspar.df}{Simulated sensitivity parameters based on observed data.}
#'   \item{res.RV}{Robustness values at specified or default representative times.}
#'   \item{summary.tables}{User-facing summary tables with clearer column names.}
#'   \item{var_names}{Confounder names used by interpretation helpers.}
#'   \item{time.info}{Time grids used for analysis and nuisance estimation.}
#' }
#'
#' @examples
#' # Simulate toy data and fit sensitivity bounds (see package vignettes for full examples)
#' \dontrun{
#' dat <- sim.data.surv.wu(n = 500, seed = 123)
#' out <- npsa_surv(
#'   time = dat$Y,
#'   event = dat$D,
#'   treat = dat$A,
#'   confounders = dat$W,
#'   rmst = FALSE
#' )
#' }
#'
#' @export
npsa_surv <- function(time, event, treat, confounders, fit.times = NULL,
                      nuisance.options = list(),
                      target.options = list(),
                      bound.options = list(),
                      rv.options = list(),
                      rmst = FALSE,
                      rmst.options = list(),
                      sens.options = list(),
                      result = NULL,
                      var_names = NULL,
                      verbose = FALSE,
                      save = FALSE) {

    # Update control parameters
    target.options <- do.call(npsa_target.options, target.options)
    bound.options <- do.call(npsa_bound.options, bound.options)
    rv.options <- do.call(npsa_rv.options, rv.options)
    rmst.options <- do.call(npsa_rmst.options, rmst.options)
    sens.options <- do.call(npsa_sens.options, sens.options)

    # Extract options
    psi.type <- target.options$psi.type
    tau.type <- target.options$tau.type
    report.times <- bound.options$report.times
    transform <- bound.options$transform
    scale <- bound.options$scale
    uniform.cutpoint <- bound.options$uniform.cutpoint
    uniform.window <- bound.options$uniform.window
    rv.times <- rv.options$rv.times
    rv.uniform.cutpoint <- rv.options$uniform.cutpoint
    rv.uniform.window <- rv.options$uniform.window
    rho <- rv.options$rho
    theta <- rv.options$theta
    fit.times.rmst <- rmst.options$fit.times.rmst
    gamma.type <- rmst.options$gamma.type
    max_gap <- rmst.options$max_gap
    tol <- rmst.options$tol
    tol1 <- rmst.options$tol1
    tol2 <- rmst.options$tol2
    pct_drop <- sens.options$pct_drop
    rep <- sens.options$rep
    senspar.df <- sens.options$senspar.df
    num_drop <- sens.options$num_drop
    senspar.save.path <- sens.options$senspar.save.path
    sens.seed <- sens.options$seed
    senspar.only <- sens.options$senspar.only

    n_var <- ncol(confounders)
    if (is.null(var_names)) {
        var_names <- colnames(confounders)
        if (is.null(var_names)) var_names <- paste0("W", seq_len(n_var))
    }
    if (length(var_names) != n_var) {
        stop("`var_names` must have one name for each confounder.")
    }
    fit.times.source <- NULL
    if (is.null(fit.times) && !is.null(result$fit.times)) {
        fit.times <- result$fit.times
        fit.times.source <- "result"
    }
    time.rst <- .get.time.info(time, event, fit.times, nuisance.options,
                               verbose = verbose)
    fit.times <- time.rst$fit.times
    nuisance.options <- time.rst$nuisance.options
    time.info <- time.rst$time.info
    if (!is.null(fit.times.source)) time.info$fit.times.source <- fit.times.source

    # Nuisance Estimation
    if (verbose) cat("Start estimating nuisances:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
    if (is.null(result) || is.null(result$nuisance)) {
        result <- .get.nuisances.est(time, event, treat, confounders, fit.times,
                                     nuisance.options = nuisance.options, verbose = verbose)

        if (save) save(result, file = "dev/result.RData")
    }

    # Observed Components Estimation
    if (verbose) cat("Start estimating target:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
    if (is.null(result$obs.comps.df)) {
        result <- .get.obs.comps(time, event, treat, result,
                                 psi.type = psi.type, tau.type = tau.type,
                                 verbose = verbose)
        if (save) save(result, file = "dev/result.RData")
    }
    if (!is.null(result$obs.comps.df$gamma)) result$obs.comps.df$gamma <- NULL
    time.info$fit.times <- result$fit.times
    time.info$eval.times <- result$nuisance$eval.times
    time.info$n.fit.times <- length(result$fit.times)
    time.info$n.eval.times <- length(result$nuisance$eval.times)

    # RMST Estimation if requested
    if (rmst) {
        fit.times.rmst.source <- "user"
        if (is.null(fit.times.rmst) && !is.null(result$fit.times.rmst)) {
            fit.times.rmst <- result$fit.times.rmst
            fit.times.rmst.source <- "result"
        } else if (is.null(fit.times.rmst)) {
            fit.times.rmst <- max(result$fit.times)
            fit.times.rmst.source <- "default"
            message("`rmst.options$fit.times.rmst` was not supplied. Using max(fit.times) = ",
                    signif(fit.times.rmst, 4),
                    " as the RMST horizon. Please check whether this RMST horizon is meaningful for your analysis, because RMST requires analyst interpretation.")
        }
        if (verbose) cat("Start estimating RMST:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
        if (is.null(result$rmst.obs)) {
            eval.times.rmst <- result$fit.times
            # cat(eval.times.rmst, "\n")
            result <- .get.rmst.obs.comps(time, event, result, fit.times.rmst, eval.times.rmst,
                                          max_gap, tol, tol1, tol2,
                                          gamma.type, verbose = verbose)
            if (save) save(result, file = "dev/result.RData")
        }
        time.info$fit.times.rmst <- result$fit.times.rmst
        time.info$fit.times.rmst.source <- fit.times.rmst.source
        rmst.options$fit.times.rmst <- result$fit.times.rmst
    }

    # Simulate sensitivity parameters if needed
    senspar.saved <- FALSE
    if (is.null(senspar.df)) {
        if (verbose) cat("Start simulating sensitivity parameters:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
        senspar.df <- .simulate.senspar(time, event, treat, confounders,
                                        fit.times = result$fit.times,
                                        psi = result$obs.comps.df$psi,
                                        tau = result$tau,
                                        S.hat.obs = result$nuisance$event.pred,
                                        g.hat.obs = result$nuisance$prop.pred,
                                        num_drop = num_drop,
                                        pct_drop = pct_drop,
                                        rep = rep,
                                        seed = sens.seed,
                                        rmst = rmst,
                                        fit.times.rmst = if (rmst) result$fit.times.rmst else NULL,
                                        gamma = if (rmst) result$gamma.est else NULL,
                                        max_gap = if (rmst) max_gap else NULL,
                                        tol = if (rmst) tol else NULL,
                                        var_names = var_names,
                                        verbose = verbose)
        if (save || !is.null(senspar.save.path)) {
            if (is.null(senspar.save.path)) senspar.save.path <- "dev/senspar.df.RData"
            dir.create(dirname(senspar.save.path), recursive = TRUE, showWarnings = FALSE)
            save(senspar.df, file = senspar.save.path)
            senspar.saved <- TRUE
        }
    } else {
        if (verbose) cat("Using user-provided sensitivity parameters:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
        senspar.df <- senspar.df
    }

    if (senspar.only) {
        if ((save || !is.null(senspar.save.path)) && !senspar.saved) {
            if (is.null(senspar.save.path)) senspar.save.path <- "dev/senspar.df.RData"
            dir.create(dirname(senspar.save.path), recursive = TRUE, showWarnings = FALSE)
            save(senspar.df, file = senspar.save.path)
        }
        out <- list(result = result,
                    senspar.df = senspar.df,
                    var_names = var_names,
                    time.info = time.info)
        out$summary.tables <- .npsa_surv_summary_tables(out)
        class(out) <- "npsa_surv"
        if (verbose) cat("Finished after sensitivity parameter simulation:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
        return(out)
    }

    # Report times and uniform window
    if (is.null(report.times)) {
        report.times <- .get.report.times(report.times, result$fit.times, default.all = TRUE,
                                        label = "report.times")
        if (is.null(rv.times)) {
            rv.times <- result$fit.times
            if (length(rv.times) > 5) {
                idx <- unique(round(seq(1, length(rv.times), length.out = 5)))
                rv.times <- rv.times[idx]
            }
        }
    } else {
        report.times <- .get.report.times(report.times, result$fit.times, default.all = TRUE,
                                        label = "report.times")
        if (is.null(rv.times)) {
            rv.times <- report.times
        }
    }
    if (!is.null(rv.times)) {
        rv.times <- .get.report.times(rv.times, result$fit.times, default.all = FALSE,
                                      label = "rv.times")
    }
    if (is.null(uniform.window)) {
        surv.0 <- .np_get_surv_object(result, trt = 0, isotonize = TRUE)
        surv.1 <- .np_get_surv_object(result, trt = 1, isotonize = TRUE)
        uniform.window <- .np_contrast_band_endpts(time[event == 1],
                                                   surv.0$surv.iso, surv.1$surv.iso,
                                                   result$fit.times, uniform.cutpoint)
        uniform.window.source <- "uniform.cutpoint"
    } else {
        uniform.window <- .np_trim_uniform_window(uniform.window, result$fit.times)
        uniform.window.source <- "uniform.window"
    }
    if (!any(result$fit.times >= uniform.window[1] & result$fit.times <= uniform.window[2])) {
        stop("No `fit.times` fall inside `uniform.window`.")
    }
    time.info$report.times <- report.times
    time.info$rv.times <- rv.times
    time.info$uniform.window <- uniform.window
    time.info$uniform.window.source <- uniform.window.source
    time.info$uniform.cutpoint <- uniform.cutpoint
    bound.times <- result$fit.times
    urv.window <- uniform.window
    urv.window.source <- uniform.window.source
    if (!is.null(rv.uniform.window)) {
        urv.window <- .np_trim_uniform_window(rv.uniform.window, result$fit.times)
        urv.window.source <- "rv.options$uniform.window"
    } else if (!is.null(rv.uniform.cutpoint)) {
        surv.0 <- .np_get_surv_object(result, trt = 0, isotonize = TRUE)
        surv.1 <- .np_get_surv_object(result, trt = 1, isotonize = TRUE)
        urv.window <- .np_contrast_band_endpts(time[event == 1],
                                               surv.0$surv.iso, surv.1$surv.iso,
                                               result$fit.times, rv.uniform.cutpoint)
        urv.window.source <- "rv.options$uniform.cutpoint"
    }

    # Observed bounds
    if (verbose) cat("Start computing observed bounds:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
    bounds.df <- .report.bounds(bound.times, result, rho = rho, rmst = rmst, transform = transform,
                                scale = scale, band.end.pts = uniform.window)

    # Bounds under sensitivity
    if (verbose) cat("Start computing sensitivity bounds:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
    bounds.df.sens <- .report.bounds(bound.times, result,
                                rho = rho,
                                sens.df.mean = senspar.df$sens.df.mean,
                                num_drop = num_drop,
                                pct_drop = pct_drop,
                                n_var = n_var,
                                rmst = rmst,
                                sens.rmst.df.mean = senspar.df$sens.rmst.df.mean,
                                transform = transform, scale = scale,
                                band.end.pts = uniform.window)

    bounds.df$bounds.df <- rbind(bounds.df$bounds.df, bounds.df.sens$bounds.df)
    if (rmst) bounds.df$bounds.df.rmst <- rbind(bounds.df$bounds.df.rmst, bounds.df.sens$bounds.df.rmst)


    # Plot if requested
    # if (plot) {
    #     if (verbose) cat("Start plotting:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
    #     plot(bounds.df$bounds.df)
    # }

    out <- list(result = result, senspar.df = senspar.df, bounds.df = bounds.df,
                var_names = var_names, time.info = time.info)

    # Robustness Values computations
    if (length(rv.times) > 0) {
        if (verbose) cat("Start computing robustness values (RV):", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
        t.lower <- urv.window[1]
        t.upper <- urv.window[2]
        uniform.test.sp0 <- NULL
        if (!is.null(result$uniform.test.sp0)) {
            stored.test <- result$uniform.test.sp0
            stored.window <- as.numeric(c(stored.test$t.lower[1], stored.test$t.upper[1]))
            if (isTRUE(all.equal(stored.window, as.numeric(urv.window))) &&
                isTRUE(all.equal(as.numeric(stored.test$theta[1]), as.numeric(theta)))) {
                uniform.test.sp0 <- stored.test
            }
        }

        out$res.RV <- .report.RV(rv.times, result, rho = rho, theta = theta,
                                 transform = transform,
                                 verbose = verbose,
                                 unif = TRUE, t.lower = t.lower, t.upper = t.upper,
                                 uniform.test.sp0 = uniform.test.sp0)
        out$res.RV$uniform.window <- urv.window
        out$res.RV$uniform.window.source <- urv.window.source
    }

    out$summary.tables <- .npsa_surv_summary_tables(out)
    class(out) <- "npsa_surv"
    if (verbose) cat("Finished:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
    return(out)
}





#' No-Unobserved-Confounding Survival Analysis
#'
#' Estimate the adjusted survival results under no unobserved confounding.
#' This user-facing wrapper reuses the same nuisance and observed-component
#' pipeline used by \code{\link{npsa_surv}()}, then reports treatment-specific
#' survival curves, the survival-difference \code{sp = 0} bounds dataframe,
#' and common CFsurvival-style contrasts not defined by the sensitivity target.
#'
#' If \code{fit.times} is not supplied, the function chooses a compact analysis
#' grid from positive observed follow-up times before the largest observed event
#' time, following the CFsurvival default idea. Most users do not need to set
#' \code{nuisance.options$eval.times}; it is an internal prediction grid and is
#' chosen automatically when omitted.
#'
#' @param time Numeric vector of event or censoring times.
#' @param event Numeric vector of event indicators (1 = event, 0 = censored).
#' @param treat Numeric vector of treatment assignment indicators (1 = treated, 0 = control).
#' @param confounders Matrix or data frame of observed confounders.
#' @param fit.times Optional numeric vector of times for survival estimation.
#'   If \code{NULL}, a compact default grid is chosen from the observed
#'   follow-up times.
#' @param nuisance.options List of options for nuisance estimation.
#' @param np.options List of options from \code{\link{np_surv.options}()}.
#' @param rmst Logical; if TRUE, estimate RMST difference under no unobserved
#'   confounding.
#' @param rmst.options List of options for RMST estimation. If \code{rmst = TRUE}
#'   and \code{fit.times.rmst} is not supplied, the default RMST horizon is
#'   \code{max(fit.times)} with a message reminding users that RMST horizons
#'   need analyst interpretation.
#' @param result Optional precomputed result object containing nuisances or observed components.
#' @param var_names Character vector of confounder variable names.
#' @param verbose Logical; if TRUE, print progress messages.
#' @param save Logical; if TRUE, save intermediate result to \code{dev/result.RData}.
#'
#' @return A list of class \code{npSurv}, including \code{time.info} with the
#'   selected analysis and nuisance time grids, and \code{summary.tables} with
#'   clean user-facing summary tables. When \code{rmst = TRUE}, the returned
#'   object also includes \code{rmst.diff.df}.
#'
#' @export
np_surv <- function(time, event, treat, confounders, fit.times = NULL,
                    nuisance.options = list(),
                    np.options = list(),
                    rmst = FALSE,
                    rmst.options = list(),
                    result = NULL,
                    var_names = NULL,
                    verbose = FALSE,
                    save = FALSE) {

    # Update control parameters
    np.options <- do.call(np_surv.options, np.options)
    rmst.options <- do.call(npsa_rmst.options, rmst.options)

    # Extract options
    psi.type <- "hybrid"
    tau.type <- "hybrid"
    report.times <- np.options$report.times
    conf.band <- np.options$conf.band
    conf.level <- np.options$conf.level
    contrasts <- np.options$contrasts
    uniform.cutpoint <- np.options$uniform.cutpoint
    uniform.window <- np.options$uniform.window
    isotonize <- np.options$isotonize
    seed <- np.options$seed
    fit.times.rmst <- rmst.options$fit.times.rmst
    gamma.type <- rmst.options$gamma.type
    max_gap <- rmst.options$max_gap
    tol <- rmst.options$tol
    tol1 <- rmst.options$tol1
    tol2 <- rmst.options$tol2

    n_var <- ncol(confounders)
    if (is.null(var_names)) {
        var_names <- colnames(confounders)
        if (is.null(var_names)) var_names <- paste0("W", seq_len(n_var))
    }
    if (length(var_names) != n_var) {
        stop("`var_names` must have one name for each confounder.")
    }
    fit.times.source <- NULL
    if (is.null(fit.times) && !is.null(result$fit.times)) {
        fit.times <- result$fit.times
        fit.times.source <- "result"
    }
    time.rst <- .get.time.info(time, event, fit.times, nuisance.options,
                               verbose = verbose)
    fit.times <- time.rst$fit.times
    nuisance.options <- time.rst$nuisance.options
    time.info <- time.rst$time.info
    if (!is.null(fit.times.source)) time.info$fit.times.source <- fit.times.source
    if (is.null(seed)) seed <- sample(1:1e8, 1)
    set.seed(seed)
    np.options$seed <- seed

    # Nuisance Estimation
    if (verbose) cat("Start estimating nuisances:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
    if (is.null(result) || is.null(result$nuisance)) {
        result <- .get.nuisances.est(time, event, treat, confounders, fit.times,
                                     nuisance.options = nuisance.options, verbose = verbose)
        if (save) save(result, file = "dev/result.RData")
    }

    # Observed Components Estimation
    if (verbose) cat("Start estimating no-unobserved-confounding survival:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
    if (is.null(result$obs.comps.df)) {
        result <- .get.obs.comps(time, event, treat, result,
                                 psi.type = psi.type, tau.type = tau.type,
                                 verbose = verbose)
        if (save) save(result, file = "dev/result.RData")
    }
    if (!is.null(result$obs.comps.df$gamma)) result$obs.comps.df$gamma <- NULL
    time.info$fit.times <- result$fit.times
    time.info$eval.times <- result$nuisance$eval.times
    time.info$n.fit.times <- length(result$fit.times)
    time.info$n.eval.times <- length(result$nuisance$eval.times)

    # RMST Estimation if requested
    if (rmst) {
        fit.times.rmst.source <- "user"
        if (is.null(fit.times.rmst) && !is.null(result$fit.times.rmst)) {
            fit.times.rmst <- result$fit.times.rmst
            fit.times.rmst.source <- "result"
        } else if (is.null(fit.times.rmst)) {
            fit.times.rmst <- max(result$fit.times)
            fit.times.rmst.source <- "default"
            message("`rmst.options$fit.times.rmst` was not supplied. Using max(fit.times) = ",
                    signif(fit.times.rmst, 4),
                    " as the RMST horizon. Please check whether this RMST horizon is meaningful for your analysis, because RMST requires analyst interpretation.")
        }
        if (verbose) cat("Start estimating RMST:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
        if (is.null(result$rmst.obs)) {
            eval.times.rmst <- result$fit.times
            result <- .get.rmst.obs.comps(time, event, result, fit.times.rmst, eval.times.rmst,
                                          max_gap, tol, tol1, tol2,
                                          gamma.type, verbose = verbose)
            if (save) save(result, file = "dev/result.RData")
        }
        time.info$fit.times.rmst <- result$fit.times.rmst
        time.info$fit.times.rmst.source <- fit.times.rmst.source
        rmst.options$fit.times.rmst <- result$fit.times.rmst
    }

    report.times <- .get.report.times(report.times, result$fit.times, default.all = TRUE,
                                    label = "report.times")
    time.info$report.times <- report.times

    # Treatment-specific survival and survival contrasts
    cf.contrasts <- setdiff(contrasts, "surv.diff")
    cf.out <- .np_report_cf_surv(time, event, treat, result,
                                 conf.band = conf.band,
                                 conf.level = conf.level,
                                 contrasts = cf.contrasts,
                                 uniform.cutpoint = uniform.cutpoint,
                                 uniform.window = uniform.window,
                                 isotonize = isotonize)

    # Survival difference is the zero-sensitivity special case of SurvNPSA.
    surv.diff.out <- .report.bounds(result$fit.times, result,
                                    rmst = FALSE,
                                    transform = TRUE,
                                    scale = TRUE,
                                    band.end.pts = cf.out$band.end.pts,
                                    conf.level = conf.level)
    cf.out$surv.diff.df <- surv.diff.out$bounds.df
    rmst.diff.df <- NULL
    rmst.summary <- NULL
    if (rmst) {
        rmst.out <- .report.bounds(result$fit.times, result,
                                   rmst = TRUE,
                                   transform = TRUE,
                                   scale = TRUE,
                                   band.end.pts = cf.out$band.end.pts,
                                   conf.level = conf.level)
        rmst.diff.df <- rmst.out$bounds.df.rmst
        rmst.summary <- .np_ci_summary(rmst.diff.df, result$fit.times.rmst)
        names(rmst.summary)[names(rmst.summary) == "surv.diff"] <- "rmst.diff"
    }

    # Uniform test for no observed survival difference
    uniform.test <- .np_uniform_test(result, time, event,
                                    uniform.cutpoint = uniform.cutpoint,
                                    uniform.window = uniform.window,
                                    conf.level = conf.level,
                                    seed = seed)
    result$uniform.test.sp0 <- uniform.test
    time.info$uniform.window <- cf.out$band.end.pts
    time.info$uniform.window.source <- if (is.null(uniform.window)) "uniform.cutpoint" else "uniform.window"
    time.info$uniform.cutpoint <- uniform.cutpoint

    ci.summary <- .np_ci_summary(cf.out$surv.diff.df, report.times)

    out <- c(list(result = result,
                  uniform.test = uniform.test,
                  ci.summary = ci.summary,
                  rmst.summary = rmst.summary,
                  report.times = report.times,
                  var_names = var_names,
                  time.info = time.info,
                  options = list(np.options = np.options,
                                 rmst.options = rmst.options)),
             cf.out)
    if (rmst) out$rmst.diff.df <- rmst.diff.df
    out$summary.tables <- .np_surv_summary_tables(out)

    class(out) <- "npSurv"
    if (verbose) cat("Finished:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
    return(out)
}

#' Options for \code{np_surv()}
#'
#' @param report.times Optional numeric vector of times to summarize and report.
#'   If \code{NULL}, all fitted times are used.
#' @param conf.band Logical; if TRUE, compute uniform confidence bands.
#' @param conf.level Desired confidence level.
#' @param contrasts Character vector of contrasts to report. Options are
#'   \code{"surv.diff"}, \code{"surv.ratio"}, \code{"risk.ratio"}, and \code{"nnt"}.
#' @param uniform.cutpoint Two probabilities used for uniform procedures. The
#'   first gives the lower event-time quantile and the second gives the upper
#'   survival threshold through \code{1 - p}, following the CFsurvival uniform
#'   band window.
#' @param uniform.window Optional exact two-time window for uniform procedures.
#'   If supplied, it is used instead of \code{uniform.cutpoint}.
#' @param isotonize Logical; if TRUE, apply CFsurvival-style isotonization to
#'   treatment-specific survival curve display and treatment-specific survival
#'   uniform bands. Pointwise confidence intervals and survival contrasts remain
#'   based on the raw one-step survival estimates and influence functions.
#' @param seed Optional integer seed for reproducible uniform bands and uniform
#'   test p-values. Use \code{NULL} to randomly choose and store a seed.
#'
#' @return A named list of options.
#'
#' @export
np_surv.options <- function(report.times = NULL, conf.band = TRUE, conf.level = 0.95,
                            contrasts = c("surv.diff", "surv.ratio", "risk.ratio", "nnt"),
                            uniform.cutpoint = c(0.01, 0.99), uniform.window = NULL,
                            isotonize = TRUE, seed = NULL) {
    if (!is.null(report.times) &&
        (!is.numeric(report.times) || any(!is.finite(report.times)) || any(report.times < 0))) {
        stop("`report.times` must be NULL or a non-negative numeric vector.")
    }
    if (length(conf.band) != 1 || !is.logical(conf.band) || is.na(conf.band)) {
        stop("`conf.band` must be TRUE or FALSE.")
    }
    if (length(conf.level) != 1 || !is.numeric(conf.level) ||
        !is.finite(conf.level) || conf.level <= 0 || conf.level >= 1) {
        stop("`conf.level` must be a number between 0 and 1.")
    }
    if (length(uniform.cutpoint) != 2 || !is.numeric(uniform.cutpoint) ||
        any(!is.finite(uniform.cutpoint)) ||
        any(uniform.cutpoint <= 0 | uniform.cutpoint >= 1) ||
        uniform.cutpoint[1] >= uniform.cutpoint[2]) {
        stop("`uniform.cutpoint` must contain two increasing numbers between 0 and 1.")
    }
    if (!is.null(uniform.window) &&
        (length(uniform.window) != 2 || !is.numeric(uniform.window) ||
         any(!is.finite(uniform.window)) || uniform.window[1] >= uniform.window[2])) {
        stop("`uniform.window` must be NULL or two increasing finite numbers.")
    }
    if (length(isotonize) != 1 || !is.logical(isotonize) || is.na(isotonize)) {
        stop("`isotonize` must be TRUE or FALSE.")
    }
    if (!is.null(seed)) {
        if (!is.numeric(seed) || length(seed) != 1 || !is.finite(seed) ||
            seed < 0 || seed > .Machine$integer.max || seed != floor(seed)) {
            stop("`seed` must be NULL or a single non-negative integer.")
        }
        seed <- as.integer(seed)
    }

    allowed <- c("surv.diff", "surv.ratio", "risk.ratio", "nnt")
    if (is.null(contrasts)) {
        contrasts <- character(0)
    } else {
        contrasts <- unique(tolower(contrasts))
        invalid <- setdiff(contrasts, allowed)
        if (length(invalid) > 0) {
            stop("Invalid contrast(s): ", paste(invalid, collapse = ", "), ".")
        }
    }
    contrasts <- unique(c("surv.diff", contrasts))

    list(report.times = report.times,
         conf.band = conf.band,
         conf.level = conf.level,
         contrasts = contrasts,
         uniform.cutpoint = uniform.cutpoint,
         uniform.window = uniform.window,
         isotonize = isotonize,
         seed = seed)
}

npsa_target.options <- function(psi.type = "hybrid", tau.type = "hybrid") {
    list(psi.type = psi.type, tau.type = tau.type)
}

npsa_bound.options <- function(report.times = NULL, transform = TRUE, scale = TRUE,
                               uniform.cutpoint = c(0.01, 0.99),
                               uniform.window = NULL) {
    if (length(uniform.cutpoint) != 2 || !is.numeric(uniform.cutpoint) ||
        any(!is.finite(uniform.cutpoint)) ||
        any(uniform.cutpoint <= 0 | uniform.cutpoint >= 1) ||
        uniform.cutpoint[1] >= uniform.cutpoint[2]) {
        stop("`uniform.cutpoint` must contain two increasing numbers between 0 and 1.")
    }
    if (!is.null(uniform.window) &&
        (length(uniform.window) != 2 || !is.numeric(uniform.window) ||
         any(!is.finite(uniform.window)) || uniform.window[1] >= uniform.window[2])) {
        stop("`uniform.window` must be NULL or two increasing finite numbers.")
    }
    list(report.times = report.times, transform = transform, scale = scale,
         uniform.cutpoint = uniform.cutpoint, uniform.window = uniform.window)
}

npsa_rv.options <- function(rv.times = NULL, uniform.cutpoint = NULL,
                            uniform.window = NULL,
                            rho = 1, theta = 0) {
    if (!is.null(uniform.cutpoint) &&
        (length(uniform.cutpoint) != 2 || !is.numeric(uniform.cutpoint) ||
         any(!is.finite(uniform.cutpoint)) ||
         any(uniform.cutpoint <= 0 | uniform.cutpoint >= 1) ||
         uniform.cutpoint[1] >= uniform.cutpoint[2])) {
        stop("`uniform.cutpoint` must contain two increasing numbers between 0 and 1.")
    }
    if (!is.null(uniform.window) &&
        (length(uniform.window) != 2 || !is.numeric(uniform.window) ||
         any(!is.finite(uniform.window)) || uniform.window[1] >= uniform.window[2])) {
        stop("`uniform.window` must be NULL or two increasing finite numbers.")
    }
    list(rv.times = rv.times, uniform.cutpoint = uniform.cutpoint,
         uniform.window = uniform.window,
         rho = rho, theta = theta)
}

npsa_rmst.options <- function(fit.times.rmst = NULL, gamma.type = "hybrid",
                              max_gap = 0.2, tol = 0.01, tol1 = 0.01, tol2 = 0.01) {
    if (!is.null(fit.times.rmst) &&
        (!is.numeric(fit.times.rmst) || any(!is.finite(fit.times.rmst)) ||
         any(fit.times.rmst <= 0))) {
        stop("'fit.times.rmst' must be NULL or positive finite values.")
    }
    list(fit.times.rmst = fit.times.rmst, gamma.type = gamma.type,
         max_gap = max_gap, tol = tol, tol1 = tol1, tol2 = tol2)
}

npsa_sens.options <- function(pct_drop = c(0.3, 0.7), rep = 10,
                              senspar.df = NULL, num_drop = NULL,
                              senspar.save.path = NULL,
                              seed = 6741,
                              senspar.only = FALSE) {
    # Use either pct_drop or num_drop. The simulation also adds d = 1 and
    # d = ceiling(0.5 * n_var) as benchmark drop sizes when possible.
    # If num_drop is supplied by itself, use exact drop sizes instead of the default pct_drop.
    if (!is.null(num_drop) && missing(pct_drop)) pct_drop <- NULL

    if (!is.null(senspar.save.path) &&
        (!is.character(senspar.save.path) || length(senspar.save.path) != 1 ||
         is.na(senspar.save.path) || !nzchar(senspar.save.path))) {
        stop("'senspar.save.path' must be a non-empty character string or NULL.")
    }
    if (length(seed) != 1 || !is.numeric(seed) || is.na(seed)) {
        stop("'seed' must be one number.")
    }
    if (length(senspar.only) != 1 || !is.logical(senspar.only) || is.na(senspar.only)) {
        stop("'senspar.only' must be TRUE or FALSE.")
    }

    list(pct_drop = pct_drop, rep = rep, senspar.df = senspar.df,
         num_drop = num_drop, senspar.save.path = senspar.save.path,
         seed = seed, senspar.only = senspar.only)
}
