#' Nonparametric Sensitivity Analysis for Survival Outcomes
#'
#' This function estimates the treatment effect on survival outcomes with
#' nonparametric sensitivity analysis under unmeasured confounding, optionally
#' including restricted mean survival time (RMST) analysis and robustness value computation.
#'
#' @param time Numeric vector of event or censoring times.
#' @param event Numeric vector of event indicators (1 = event, 0 = censored).
#' @param treat Numeric vector of treatment assignment indicators (1 = treated, 0 = control).
#' @param confounders Matrix or data frame of observed confounders (observed covariates).
#' @param fit.times Numeric vector of times at which nuisance estimators are fit.
#' @param nuisance.options List of options for nuisance estimation.
#' @param target.options List of options for target parameter estimation.
#' @param bound.options List of options for reporting pointwise and uniform bounds.
#' @param rv.options List of options for robustness value computation.
#' @param rmst Logical; if TRUE, estimate RMST and its bounds inference as well.
#' @param rmst.options List of options for RMST estimation.
#' @param sens.options List of options for sensitivity parameter simulation.
#' @param result Optional precomputed result object (e.g., containing nuisances).
#' @param var_names Character vector of confounder variable names.
#' @param plot Logical; if TRUE, automatically plot bounds after estimation.
#' @param verbose Logical; if TRUE, print system timestamps for each estimation step.
#'
#' @return A list of class \code{npsa_surv} containing:
#' \describe{
#'   \item{result}{Estimated observable components and IFs, such as observed survival differences and rmst differences.}
#'   \item{bounds.df}{Estimated bounds on survival contrasts over time.}
#'   \item{senspar.df}{Simulated sensitivity parameters based on observed data.}
#'   \item{res.RV}{Robustness values at specified times (if \code{rv.times} is given).}
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
#'   fit.times = seq(0.5, 10, by = 0.5)
#' )
#' }
#'
#' @export
npsa_surv <- function(time, event, treat, confounders, fit.times,
                      nuisance.options = list(),
                      target.options = list(),
                      bound.options = list(),
                      rv.options = list(),
                      rmst = TRUE,
                      rmst.options = list(),
                      sens.options = list(),
                      result = NULL,
                      var_names = NULL,
                      plot = TRUE,
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
    plot.times <- bound.options$plot.times
    rv.times <- rv.options$rv.times
    uniform.cutpoint <- rv.options$uniform.cutpoint
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

    n_var <- ncol(confounders)

    ## --- Estimation Process ---

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
        result <- .get.obs.comps(time, event, treat, result, psi.type = psi.type, verbose = verbose)
        if (save) save(result, file = "dev/result.RData")
    }

    # RMST Estimation if requested
    if (rmst) {
        if (is.null(fit.times.rmst)) stop("Must specify 'fit.times.rmst' when rmst = TRUE.")
        if (verbose) cat("Start estimating RMST:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
        if (is.null(result$rmst.obs)) {
            eval.times.rmst <- result$fit.times
            cat(eval.times.rmst, "\n")
            result <- .get.rmst.obs.comps(time, event, result, fit.times.rmst, eval.times.rmst,
                                          max_gap, tol, tol1, tol2,
                                          gamma.type, verbose = verbose)
            if (save) save(result, file = "dev/result.RData")
        }
    }

    # Observed bounds
    if (verbose) cat("Start computing observed bounds:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
    bounds.df <- .report.bounds(plot.times, result, rmst = rmst)

    # Simulate sensitivity parameters if needed
    if (verbose) cat("Start simulating sensitivity parameters:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
    if (is.null(senspar.df)) {
        senspar.df <- .simulate.senspar(time, event, treat, confounders,
                                        fit.times = result$fit.times,
                                        psi = result$obs.comps.df$psi,
                                        tau = result$tau,
                                        S.hat.obs = result$nuisance$event.pred,
                                        g.hat.obs = result$nuisance$prop.pred,
                                        pct_drop = pct_drop,
                                        rep = rep,
                                        rmst = rmst,
                                        fit.times.rmst = if (rmst) result$fit.times.rmst else NULL,
                                        gamma = if (rmst) result$gamma.est else NULL,
                                        max_gap = if (rmst) max_gap else NULL,
                                        tol = if (rmst) tol else NULL)
        if (save) save(senspar.df, file = "dev/senspar.df.RData")
    } else {
        senspar.df <- senspar.df
    }

    # Bounds under sensitivity
    if (verbose) cat("Start computing sensitivity bounds:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
    bounds.df.sens <- .report.bounds(plot.times, result,
                                sens.df.mean = senspar.df$sens.df.mean,
                                num_drop = num_drop,
                                pct_drop = pct_drop,
                                n_var = n_var,
                                rmst = rmst,
                                sens.rmst.df.mean = senspar.df$sens.rmst.df.mean)

    bounds.df$bounds.df <- rbind(bounds.df$bounds.df, bounds.df.sens$bounds.df)
    if (rmst) bounds.df$bounds.df.rmst <- rbind(bounds.df$bounds.df.rmst, bounds.df.sens$bounds.df.rmst)


    # Plot if requested
    if (plot) {
        if (verbose) cat("Start plotting:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
        plot(bounds.df$bounds.df)
    }

    out <- list(result = result, senspar.df = senspar.df, bounds.df = bounds.df)

    # Robustness Values computations
    if (!is.null(rv.times)) {
        if (verbose) cat("Start computing robustness values (RV):", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
        q.01 <- quantile(time[event == 1], uniform.cutpoint[1])
        q.99 <- quantile(time[event == 1], uniform.cutpoint[2])
        cat(q.01, q.99, "\n")

        out$res.RV <- .report.RV(rv.times, result, unif = TRUE, q.01 = q.01, q.99 = q.99)
    }

    class(out) <- "npsa_surv"
    if (verbose) cat("Finished:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
    return(out)
}



npsa_target.options <- function(psi.type = "hybrid") {
    list(psi.type = psi.type)
}

npsa_bound.options <- function(plot.times = c(0.5, 0.8, 1.2)) {
    list(plot.times = plot.times)
}

npsa_rv.options <- function(rv.times = c(0.3, 0.5, 1.2), uniform.cutpoint = c(0.01, 0.99)) {
    list(rv.times = rv.times, uniform.cutpoint = uniform.cutpoint)
}

npsa_rmst.options <- function(fit.times.rmst = c(0.5, 0.7, 2), gamma.type = "hybrid",
                              max_gap = 0.2, tol = 0.01, tol1 = 0.01, tol2 = 0.01) {
    list(fit.times.rmst = fit.times.rmst, gamma.type = gamma.type,
         max_gap = max_gap, tol = tol, tol1 = tol1, tol2 = tol2)
}

npsa_sens.options <- function(pct_drop = c(0.3, 0.7), rep = 10,
                              senspar.df = NULL, num_drop = NULL) {
    list(pct_drop = pct_drop, rep = rep, senspar.df = senspar.df, num_drop = num_drop)
}





