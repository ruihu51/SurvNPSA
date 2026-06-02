if (getRversion() >= "2.15.1") {
    utils::globalVariables(c(".data", "time", "surv", "trt", "times",
                             "theta.obs", "effect.lower", "effect.upper",
                             "ptwise.lower", "ptwise.upper", "unif.lower",
                             "unif.upper", "uniform.lower", "uniform.upper"))
}

#' No-Unobserved-Confounding Survival Analysis
#'
#' Estimate the adjusted survival results under no unobserved confounding.
#' This user-facing wrapper reuses the same nuisance and observed-component
#' pipeline used by \code{\link{npsa_surv}()}, then reports treatment-specific
#' survival curves and common CFsurvival-style contrasts.
#'
#' @param time Numeric vector of event or censoring times.
#' @param event Numeric vector of event indicators (1 = event, 0 = censored).
#' @param treat Numeric vector of treatment assignment indicators (1 = treated, 0 = control).
#' @param confounders Matrix or data frame of observed confounders.
#' @param fit.times Numeric vector of times for survival estimation.
#' @param nuisance.options List of options for nuisance estimation.
#' @param target.options List of options for target parameter estimation.
#' @param np.options List of options from \code{\link{np_surv.options}()}.
#' @param rmst Logical; if TRUE, estimate RMST using the existing SurvNPSA RMST pipeline.
#' @param rmst.options List of options for RMST estimation.
#' @param result Optional precomputed result object containing nuisances or observed components.
#' @param var_names Character vector of confounder variable names.
#' @param verbose Logical; if TRUE, print progress messages.
#' @param save Logical; if TRUE, save intermediate result to \code{dev/result.RData}.
#'
#' @return A list of class \code{npSurv}.
#'
#' @export
np_surv <- function(time, event, treat, confounders, fit.times,
                    nuisance.options = list(),
                    target.options = list(),
                    np.options = list(),
                    rmst = FALSE,
                    rmst.options = list(),
                    result = NULL,
                    var_names = NULL,
                    verbose = FALSE,
                    save = FALSE) {

    target.options <- do.call(npsa_target.options, target.options)
    np.options <- do.call(np_surv.options, np.options)
    rmst.options <- do.call(npsa_rmst.options, rmst.options)

    psi.type <- target.options$psi.type
    plot.times <- np.options$plot.times
    conf.band <- np.options$conf.band
    conf.level <- np.options$conf.level
    contrasts <- np.options$contrasts
    uniform.cutpoint <- np.options$uniform.cutpoint
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
    if (!is.null(seed)) {
        has.seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
        if (has.seed) old.seed <- get(".Random.seed", envir = .GlobalEnv)
        on.exit({
            if (has.seed) {
                assign(".Random.seed", old.seed, envir = .GlobalEnv)
            } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
                rm(".Random.seed", envir = .GlobalEnv)
            }
        }, add = TRUE)
        set.seed(seed)
    }

    if (verbose) cat("Start estimating nuisances:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
    if (is.null(result) || is.null(result$nuisance)) {
        result <- .get.nuisances.est(time, event, treat, confounders, fit.times,
                                     nuisance.options = nuisance.options, verbose = verbose)
        if (save) save(result, file = "dev/result.RData")
    }

    if (verbose) cat("Start estimating no-unobserved-confounding survival:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
    if (is.null(result$obs.comps.df)) {
        result <- .get.obs.comps(time, event, treat, result, psi.type = psi.type, verbose = verbose)
        if (save) save(result, file = "dev/result.RData")
    }

    if (rmst) {
        if (is.null(fit.times.rmst)) stop("Must specify 'fit.times.rmst' when rmst = TRUE.")
        if (verbose) cat("Start estimating RMST:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
        if (is.null(result$rmst.obs)) {
            eval.times.rmst <- result$fit.times
            result <- .get.rmst.obs.comps(time, event, result, fit.times.rmst, eval.times.rmst,
                                          max_gap, tol, tol1, tol2,
                                          gamma.type, verbose = verbose)
            if (save) save(result, file = "dev/result.RData")
        }
    }

    if (is.null(plot.times)) plot.times <- result$fit.times
    plot.times <- plot.times[plot.times >= min(result$fit.times) &
                                 plot.times <= max(result$fit.times)]
    if (length(plot.times) == 0) stop("No `plot.times` remain within the fitted time range.")

    cf.out <- .np_report_cf_surv(time, event, treat, result,
                                 conf.band = conf.band,
                                 conf.level = conf.level,
                                 contrasts = contrasts,
                                 uniform.cutpoint = uniform.cutpoint,
                                 isotonize = isotonize)

    bounds.df <- .report.bounds(plot.times, result, rmst = rmst,
                                transform = TRUE, scale = TRUE)

    uniform.test <- .np_uniform_test(result, time, event,
                                    uniform.cutpoint = uniform.cutpoint,
                                    conf.level = conf.level)

    ci.summary <- .np_ci_summary(cf.out$surv.diff.df, plot.times)

    out <- c(list(result = result,
                  bounds.df = bounds.df,
                  uniform.test = uniform.test,
                  ci.summary = ci.summary,
                  plot.times = plot.times,
                  var_names = var_names,
                  options = list(target.options = target.options,
                                 np.options = np.options,
                                 rmst = rmst,
                                 rmst.options = rmst.options)),
             cf.out)

    class(out) <- "npSurv"
    if (verbose) cat("Finished:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
    return(out)
}

#' Options for \code{np_surv()}
#'
#' @param plot.times Optional numeric vector of times to summarize and plot.
#' @param conf.band Logical; if TRUE, compute uniform confidence bands.
#' @param conf.level Desired confidence level.
#' @param contrasts Character vector of contrasts to report. Options are
#'   \code{"surv.diff"}, \code{"surv.ratio"}, \code{"risk.ratio"}, and \code{"nnt"}.
#' @param uniform.cutpoint Two probabilities used to choose the time window for
#'   the no-confounding uniform test and CF-style uniform bands.
#' @param isotonize Logical; if TRUE, apply isotonization to treatment-specific
#'   survival curves and survival bands.
#' @param seed Optional integer seed for reproducible uniform bands and uniform
#'   test p-values. Use \code{NULL} to leave the random seed unchanged.
#'
#' @return A named list of options.
#'
#' @export
np_surv.options <- function(plot.times = NULL, conf.band = TRUE, conf.level = 0.95,
                            contrasts = c("surv.diff", "surv.ratio", "risk.ratio", "nnt"),
                            uniform.cutpoint = c(0.01, 0.99), isotonize = TRUE,
                            seed = NULL) {
    if (!is.null(plot.times) &&
        (!is.numeric(plot.times) || any(!is.finite(plot.times)) || any(plot.times < 0))) {
        stop("`plot.times` must be NULL or a non-negative numeric vector.")
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

    list(plot.times = plot.times,
         conf.band = conf.band,
         conf.level = conf.level,
         contrasts = contrasts,
         uniform.cutpoint = uniform.cutpoint,
         isotonize = isotonize,
         seed = seed)
}

.np_report_cf_surv <- function(time, event, treat, result, conf.band=TRUE,
                               conf.level=.95, contrasts=c("surv.diff", "surv.ratio"),
                               uniform.cutpoint=c(0.01, 0.99), isotonize=TRUE) {
    fit.times <- result$fit.times
    surv.0 <- .np_get_surv_object(result, trt = 0, isotonize = isotonize)
    surv.1 <- .np_get_surv_object(result, trt = 1, isotonize = isotonize)

    band.end.pts.0 <- .np_band_endpts(time[event == 1 & treat == 0], uniform.cutpoint, fit.times)
    band.end.pts.1 <- .np_band_endpts(time[event == 1 & treat == 1], uniform.cutpoint, fit.times)
    band.end.pts <- .np_band_endpts(time[event == 1], uniform.cutpoint, fit.times)

    surv.df.0 <- .np_surv_df_one(fit.times, surv.0, trt = 0, conf.band = conf.band,
                                 band.end.pts = band.end.pts.0,
                                 conf.level = conf.level,
                                 isotonize = isotonize)
    surv.df.1 <- .np_surv_df_one(fit.times, surv.1, trt = 1, conf.band = conf.band,
                                 band.end.pts = band.end.pts.1,
                                 conf.level = conf.level,
                                 isotonize = isotonize)

    out <- list(surv.df = rbind(surv.df.0$surv.df, surv.df.1$surv.df),
                surv.0.unif.ew.quant = surv.df.0$unif.ew.quant,
                surv.0.unif.logit.quant = surv.df.0$unif.logit.quant,
                surv.1.unif.ew.quant = surv.df.1$unif.ew.quant,
                surv.1.unif.logit.quant = surv.df.1$unif.logit.quant)

    if ("surv.diff" %in% contrasts) {
        out <- c(out, .np_surv.difference(fit.times = fit.times,
                                          surv.0 = surv.0$surv,
                                          surv.1 = surv.1$surv,
                                          IF.vals.0 = surv.0$IF.vals,
                                          IF.vals.1 = surv.1$IF.vals,
                                          conf.band = conf.band,
                                          band.end.pts = band.end.pts,
                                          conf.level = conf.level))
    }
    if ("surv.ratio" %in% contrasts) {
        out <- c(out, .np_surv.ratio(fit.times = fit.times,
                                     surv.0 = surv.0$surv,
                                     surv.1 = surv.1$surv,
                                     IF.vals.0 = surv.0$IF.vals,
                                     IF.vals.1 = surv.1$IF.vals,
                                     conf.band = conf.band,
                                     band.end.pts = band.end.pts,
                                     conf.level = conf.level))
    }
    if ("risk.ratio" %in% contrasts) {
        out <- c(out, .np_risk.ratio(fit.times = fit.times,
                                     surv.0 = surv.0$surv,
                                     surv.1 = surv.1$surv,
                                     IF.vals.0 = surv.0$IF.vals,
                                     IF.vals.1 = surv.1$IF.vals,
                                     conf.band = conf.band,
                                     band.end.pts = band.end.pts,
                                     conf.level = conf.level))
    }
    if ("nnt" %in% contrasts) {
        out <- c(out, .np_nnt(fit.times = fit.times,
                              surv.0 = surv.0$surv,
                              surv.1 = surv.1$surv,
                              IF.vals.0 = surv.0$IF.vals,
                              IF.vals.1 = surv.1$IF.vals,
                              conf.band = conf.band,
                              band.end.pts = band.end.pts,
                              conf.level = conf.level))
    }

    out$band.end.pts <- band.end.pts
    return(out)
}

.np_get_surv_object <- function(result, trt, isotonize=TRUE) {
    surv.name <- paste0("surv.", trt)
    if (!is.null(result[[surv.name]]) && !is.null(result[[surv.name]]$surv)) {
        surv <- result[[surv.name]]
    } else {
        surv <- list(times = result$fit.times,
                     surv = result$obs.comps.df[[surv.name]],
                     IF.vals = result[[paste0("IF.vals.", trt)]])
    }

    if (isotonize) {
        surv$surv.iso <- NA
        surv$surv.iso[!is.na(surv$surv)] <-
            1 - stats::isoreg(surv$times[!is.na(surv$surv)], 1-surv$surv[!is.na(surv$surv)])$yf
    } else {
        surv$surv.iso <- surv$surv
    }
    return(surv)
}

.np_surv_df_one <- function(fit.times, surv, trt, conf.band=TRUE,
                            band.end.pts=c(0, Inf), conf.level=.95, isotonize=TRUE) {
    c.int <- .np_surv.confints(fit.times, surv$surv, surv$IF.vals,
                               conf.band = conf.band,
                               band.end.pts = band.end.pts,
                               conf.level = conf.level,
                               isotonize = isotonize)
    surv.df <- data.frame(time = c(0, fit.times),
                          trt = trt,
                          surv = c(1, surv$surv.iso))
    surv.df$se <- c(0, c.int$res$se)
    surv.df$se.logit <- c(0, c.int$res$se.logit)
    surv.df$ptwise.lower <- c(1, c.int$res$ptwise.lower)
    surv.df$ptwise.upper <- c(1, c.int$res$ptwise.upper)
    surv.df$ptwise.logit.lower <- c(1, c.int$res$ptwise.logit.lower)
    surv.df$ptwise.logit.upper <- c(1, c.int$res$ptwise.logit.upper)
    surv.df$unif.ew.lower <- c(1, c.int$res$unif.ew.lower)
    surv.df$unif.ew.upper <- c(1, c.int$res$unif.ew.upper)
    surv.df$unif.logit.lower <- c(1, c.int$res$unif.logit.lower)
    surv.df$unif.logit.upper <- c(1, c.int$res$unif.logit.upper)

    return(list(surv.df = surv.df,
                unif.ew.quant = c.int$unif.ew.quant,
                unif.logit.quant = c.int$unif.logit.quant))
}

.np_band_endpts <- function(event.times, uniform.cutpoint, fit.times) {
    event.times <- event.times[is.finite(event.times)]
    if (length(event.times) == 0) return(c(min(fit.times), max(fit.times)))
    out <- as.numeric(stats::quantile(event.times, uniform.cutpoint, na.rm = TRUE))
    out[1] <- max(out[1], min(fit.times))
    out[2] <- min(out[2], max(fit.times))
    if (out[1] >= out[2]) out <- c(min(fit.times), max(fit.times))
    return(out)
}

.np_uniform_test <- function(result, time, event, uniform.cutpoint=c(0.01, 0.99),
                             conf.level=.95, theta=0) {
    band.end.pts <- .np_band_endpts(time[event == 1], uniform.cutpoint, result$fit.times)
    unif.idx <- which(result$fit.times >= band.end.pts[1] & result$fit.times <= band.end.pts[2])
    if (length(unif.idx) == 0) {
        return(data.frame(t.lower = band.end.pts[1], t.upper = band.end.pts[2],
                          n.time = 0, theta = theta, p.value = NA,
                          reject.no.effect = NA))
    }

    theta.obs <- result$obs.comps.df$theta.obs[unif.idx]
    IF.vals.theta.obs <- result$IF.vals.theta.obs[,unif.idx, drop = FALSE]
    n <- nrow(IF.vals.theta.obs)
    epsilon <- .estimate.limit.dist(IF.vals = IF.vals.theta.obs)
    test.stat <- sqrt(n) * max(abs(theta.obs - theta))
    dist.null <- apply(epsilon, 1, function(x) max(abs(x)))
    pvalue <- mean(dist.null > test.stat)

    data.frame(t.lower = band.end.pts[1],
               t.upper = band.end.pts[2],
               n.time = length(unif.idx),
               theta = theta,
               test.stat = test.stat,
               p.value = pvalue,
               reject.no.effect = pvalue < 1 - conf.level)
}

.np_ci_summary <- function(surv.diff.df, plot.times) {
    idx <- sapply(plot.times, function(x) which.min(abs(surv.diff.df$time - x)))
    idx <- unique(idx)
    out <- surv.diff.df[idx, c("time", "surv.diff", "ptwise.lower", "ptwise.upper", "ptwise.pval")]
    out$ci.includes.0 <- out$ptwise.lower <= 0 & out$ptwise.upper >= 0
    rownames(out) <- NULL
    return(out)
}

#' Summarize No-Unobserved-Confounding Survival Results
#'
#' @param object An object returned by \code{\link{np_surv}()}.
#' @param digits Number of digits for printing.
#' @param ... Additional arguments.
#'
#' @return Invisibly returns \code{object}.
#'
#' @export
#' @method summary npSurv
summary.npSurv <- function(object, digits = 3, ...) {
    cat("No-Unobserved-Confounding Survival Report\n")
    cat("-----------------------------------------\n")
    cat("\nSurvival difference summary:\n")
    tbl <- object$ci.summary
    num.cols <- sapply(tbl, is.numeric)
    tbl[, num.cols] <- lapply(tbl[, num.cols, drop = FALSE], function(x) round(x, digits))
    print(tbl, row.names = FALSE)

    if (!is.null(object$uniform.test)) {
        cat("\nUniform no-effect test:\n")
        tbl <- object$uniform.test
        num.cols <- sapply(tbl, is.numeric)
        tbl[, num.cols] <- lapply(tbl[, num.cols, drop = FALSE], function(x) round(x, digits))
        print(tbl, row.names = FALSE)
    }
    invisible(object)
}

#' Plot No-Unobserved-Confounding Survival Results
#'
#' @param x An object returned by \code{\link{np_surv}()}.
#' @param type Plot type. Options are \code{"surv"}, \code{"surv.diff"},
#'   \code{"surv.ratio"}, \code{"risk.ratio"}, \code{"nnt"}, and \code{"bounds"}.
#' @param band Band type for treatment-specific survival curves.
#' @param ... Additional arguments.
#'
#' @return A \code{ggplot} object.
#'
#' @export
#' @method plot npSurv
plot.npSurv <- function(x, type = c("surv", "surv.diff", "surv.ratio", "risk.ratio", "nnt", "bounds"),
                        band = c("logit", "equal-width", "pointwise", "none"), ...) {
    type <- match.arg(type)
    band <- match.arg(band)

    if (type == "surv") {
        return(.plot.np_surv_curves(x$surv.df, band = band))
    }
    if (type == "bounds") {
        return(.plot.np_bounds(x$bounds.df$bounds.df))
    }
    return(.plot.np_contrast(x, type = type))
}

.plot.np_surv_curves <- function(df, band = "logit") {
    lower <- upper <- NULL
    if (band == "logit") {
        lower <- "unif.logit.lower"
        upper <- "unif.logit.upper"
    } else if (band == "equal-width") {
        lower <- "unif.ew.lower"
        upper <- "unif.ew.upper"
    } else if (band == "pointwise") {
        lower <- "ptwise.lower"
        upper <- "ptwise.upper"
    }

    p <- ggplot2::ggplot(df, ggplot2::aes(x = time, y = surv, color = as.factor(trt), group = trt)) +
        ggplot2::geom_step() +
        ggplot2::scale_color_manual(values = c("0" = "#0072B2", "1" = "#D55E00"),
                                    labels = c("0" = "Control", "1" = "Treatment")) +
        ggplot2::labs(color = "Treatment") +
        ggplot2::xlab("Time") +
        ggplot2::ylab("Treatment-specific survival") +
        ggplot2::coord_cartesian(ylim = c(0, 1)) +
        ggplot2::theme_bw() +
        ggplot2::theme(legend.position = "bottom",
                       legend.title = ggplot2::element_blank(),
                       panel.grid.minor = ggplot2::element_blank())

    if (!is.null(lower) && lower %in% names(df) && upper %in% names(df)) {
        p <- p +
            ggplot2::geom_step(ggplot2::aes(y = .data[[lower]]), linetype = "dashed", na.rm = TRUE) +
            ggplot2::geom_step(ggplot2::aes(y = .data[[upper]]), linetype = "dashed", na.rm = TRUE)
    }
    return(p)
}

.plot.np_contrast <- function(x, type) {
    map <- list(
        "surv.diff" = list(df = x$surv.diff.df, est = "surv.diff",
                           ylab = "Survival difference (treatment - control)", ref = 0),
        "surv.ratio" = list(df = x$surv.ratio.df, est = "surv.ratio",
                            ylab = "Survival ratio (treatment / control)", ref = 1),
        "risk.ratio" = list(df = x$risk.ratio.df, est = "risk.ratio",
                            ylab = "Risk ratio (treatment / control)", ref = 1),
        "nnt" = list(df = x$nnt.df, est = "nnt",
                     ylab = "Number needed to treat", ref = NA)
    )
    info <- map[[type]]
    if (is.null(info$df)) stop("Requested contrast was not computed.")
    df <- info$df

    p <- ggplot2::ggplot(df, ggplot2::aes(x = time)) +
        ggplot2::geom_line(ggplot2::aes(y = .data[[info$est]], color = "Estimate"), na.rm = TRUE) +
        ggplot2::geom_line(ggplot2::aes(y = ptwise.lower, color = "Pointwise CI"), linetype = "dashed", na.rm = TRUE) +
        ggplot2::geom_line(ggplot2::aes(y = ptwise.upper, color = "Pointwise CI"), linetype = "dashed", na.rm = TRUE) +
        ggplot2::geom_line(ggplot2::aes(y = unif.lower, color = "Uniform Band"), linetype = "longdash", na.rm = TRUE) +
        ggplot2::geom_line(ggplot2::aes(y = unif.upper, color = "Uniform Band"), linetype = "longdash", na.rm = TRUE) +
        ggplot2::scale_color_manual(values = c("Estimate" = "black",
                                               "Pointwise CI" = "#0072B2",
                                               "Uniform Band" = "#009E73")) +
        ggplot2::xlab("Time") +
        ggplot2::ylab(info$ylab) +
        ggplot2::theme_bw() +
        ggplot2::theme(legend.position = "bottom",
                       legend.title = ggplot2::element_blank(),
                       panel.grid.minor = ggplot2::element_blank())

    if (is.finite(info$ref)) {
        p <- p + ggplot2::geom_hline(yintercept = info$ref, color = "grey45", linetype = "dotted")
    }
    return(p)
}

.plot.np_bounds <- function(df) {
    if ("ptwise.trans.lower" %in% names(df)) {
        df$ptwise.lower <- df$ptwise.trans.lower
        df$ptwise.upper <- df$ptwise.trans.upper
        df$uniform.lower <- df$uniform.trans.lower
        df$uniform.upper <- df$uniform.trans.upper
    } else {
        df$ptwise.lower <- df$ptwise.bounds.lower
        df$ptwise.upper <- df$ptwise.bounds.upper
        df$uniform.lower <- df$uniform.bounds.lower
        df$uniform.upper <- df$uniform.bounds.upper
    }

    ggplot2::ggplot(df, ggplot2::aes(x = times)) +
        ggplot2::geom_line(ggplot2::aes(y = theta.obs, color = "Observed Effect", linetype = "Observed Effect")) +
        ggplot2::geom_line(ggplot2::aes(y = effect.lower, color = "Lower Effect Bound", linetype = "Lower Effect Bound")) +
        ggplot2::geom_line(ggplot2::aes(y = effect.upper, color = "Upper Effect Bound", linetype = "Upper Effect Bound")) +
        ggplot2::geom_line(ggplot2::aes(y = ptwise.lower, color = "Pointwise CI", linetype = "Pointwise CI")) +
        ggplot2::geom_line(ggplot2::aes(y = ptwise.upper, color = "Pointwise CI", linetype = "Pointwise CI")) +
        ggplot2::geom_line(ggplot2::aes(y = uniform.lower, color = "Uniform Band", linetype = "Uniform Band")) +
        ggplot2::geom_line(ggplot2::aes(y = uniform.upper, color = "Uniform Band", linetype = "Uniform Band")) +
        ggplot2::scale_color_manual(values = c(
            "Observed Effect" = "black",
            "Lower Effect Bound" = "#D55E00",
            "Upper Effect Bound" = "#009E73",
            "Pointwise CI" = "#0072B2",
            "Uniform Band" = "#CC79A7"
        )) +
        ggplot2::scale_linetype_manual(values = c(
            "Observed Effect" = "solid",
            "Lower Effect Bound" = "dashed",
            "Upper Effect Bound" = "dotdash",
            "Pointwise CI" = "twodash",
            "Uniform Band" = "longdash"
        )) +
        ggplot2::xlab("Time") +
        ggplot2::ylab("Survival difference (treatment - control)") +
        ggplot2::labs(color = "Type", linetype = "Type") +
        ggplot2::theme_bw() +
        ggplot2::theme(legend.position = "bottom",
                       legend.title = ggplot2::element_blank(),
                       panel.grid.minor = ggplot2::element_blank())
}
