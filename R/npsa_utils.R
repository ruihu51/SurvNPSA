#' Check default analysis times
#'
#' This helper shows the default time grid that \code{\link{np_surv}()} and
#' \code{\link{npsa_surv}()} would use when \code{fit.times} is not provided.
#' It does not fit nuisance models or estimate treatment effects.
#'
#' @param time Numeric vector of event or censoring times.
#' @param event Numeric vector of event indicators (1 = event, 0 = censored).
#' @param fit.times Optional user-provided analysis times.
#' @param nuisance.options Optional nuisance options. If
#'   \code{nuisance.options$eval.times} is supplied, it is preserved.
#' @param plot.times Optional reporting times. If \code{NULL}, all fitted times
#'   are used.
#' @param rv.times Optional RV/MIRV times. If \code{NULL}, \code{plot.times}
#'   are used when supplied; otherwise about five representative fitted times
#'   are used.
#' @param rmst Logical; if TRUE, preview RMST horizon settings.
#' @param fit.times.rmst Optional RMST horizon times. If \code{NULL} and
#'   \code{rmst = TRUE}, the preview uses \code{max(fit.times)} and reminds
#'   users that RMST horizons need analyst interpretation.
#' @param max.fit.times Maximum number of automatic fitted times.
#' @param max.eval.times Maximum number of automatic nuisance prediction times.
#' @param G.cutoff Practical reverse-KM censoring-support cutoff for automatic
#'   fitted times.
#' @param verbose Logical; if TRUE, print a short message when automatic times
#'   are chosen.
#'
#' @return A list of class \code{npsa_times} containing \code{fit.times},
#'   \code{eval.times}, \code{plot.times}, \code{rv.times}, \code{time.info},
#'   and \code{censor.df}.
#'
#' @examples
#' \dontrun{
#' time.out <- npsa_times(time, event)
#' summary(time.out)
#' plot(time.out)
#'
#' time.out <- npsa_times(time, event,
#'                        fit.times = seq(0.1, 1.2, by = 0.1),
#'                        plot.times = c(0.2, 0.6, 1.0),
#'                        rv.times = c(0.2, 0.6, 1.0))
#' summary(time.out)
#'
#' time.out <- npsa_times(time, event,
#'                        fit.times = c(0.2, 0.4, 0.6, 0.8))
#' summary(time.out)
#' }
#'
#' @export
npsa_times <- function(time, event, fit.times = NULL, nuisance.options = list(),
                       plot.times = NULL, rv.times = NULL,
                       rmst = FALSE, fit.times.rmst = NULL,
                       max.fit.times = 50, max.eval.times = 200,
                       G.cutoff = 0.05, verbose = FALSE) {

    time.rst <- .get.time.info(time, event, fit.times, nuisance.options,
                               max.fit.times = max.fit.times,
                               max.eval.times = max.eval.times,
                               G.cutoff = G.cutoff,
                               verbose = verbose)

    fit.times <- time.rst$fit.times
    eval.times <- time.rst$nuisance.options$eval.times

    if (is.null(plot.times)) {
        plot.times <- .get.report.times(plot.times, fit.times, default.all = TRUE,
                                        label = "plot.times")
        if (is.null(rv.times)) {
            rv.times <- fit.times
            if (length(rv.times) > 5) {
                idx <- unique(round(seq(1, length(rv.times), length.out = 5)))
                rv.times <- rv.times[idx]
            }
        }
    } else {
        plot.times <- .get.report.times(plot.times, fit.times, default.all = TRUE,
                                        label = "plot.times")
        if (is.null(rv.times)) {
            rv.times <- plot.times
        }
    }
    if (!is.null(rv.times)) {
        rv.times <- .get.report.times(rv.times, fit.times, default.all = FALSE,
                                      label = "rv.times")
    }

    rmst.message <- NULL
    fit.times.rmst.source <- NULL
    if (rmst) {
        if (is.null(fit.times.rmst)) {
            fit.times.rmst <- max(fit.times)
            fit.times.rmst.source <- "default"
            rmst.message <- paste0("RMST horizon was not supplied. Using max(fit.times) = ",
                                   signif(fit.times.rmst, 4),
                                   ". Please check whether this RMST horizon is meaningful for your analysis, because RMST requires analyst interpretation.")
        } else {
            if (!is.numeric(fit.times.rmst) || any(!is.finite(fit.times.rmst)) || any(fit.times.rmst <= 0)) {
                stop("`fit.times.rmst` must contain positive finite values.")
            }
            if (any(fit.times.rmst > max(fit.times))) {
                message("Some fit.times.rmst > max(fit.times) - removed for RMST preview.")
                fit.times.rmst <- fit.times.rmst[fit.times.rmst <= max(fit.times)]
            }
            fit.times.rmst <- sort(unique(fit.times.rmst))
            if (length(fit.times.rmst) == 0) {
                stop("No `fit.times.rmst` remain within the fitted time range.")
            }
            fit.times.rmst.source <- "user"
        }
    }

    out <- list(fit.times = fit.times,
                eval.times = eval.times,
                plot.times = plot.times,
                rv.times = rv.times,
                rmst = rmst,
                fit.times.rmst = fit.times.rmst,
                time.info = time.rst$time.info,
                censor.df = time.rst$censor.df)
    out$time.info$fit.times.rmst <- fit.times.rmst
    out$time.info$fit.times.rmst.source <- fit.times.rmst.source
    out$time.info$rmst.message <- rmst.message
    class(out) <- "npsa_times"
    return(out)
}

#' Summarize default analysis times
#'
#' @param object An object returned by \code{\link{npsa_times}()}.
#' @param digits Number of digits for printing.
#' @param ... Additional arguments.
#'
#' @return Invisibly returns \code{object}.
#'
#' @export
summary.npsa_times <- function(object, digits = 4, ...) {
    .show.times <- function(x) {
        paste(signif(x, digits), collapse = ", ")
    }
    cat("Default time setting\n")
    cat("fit.times:", length(object$fit.times), "from",
        signif(min(object$fit.times), digits), "to",
        signif(max(object$fit.times), digits), "\n")
    if (length(object$fit.times) <= 10) {
        cat("fit.times values:", .show.times(object$fit.times), "\n")
    }
    cat("eval.times:", length(object$eval.times), "from",
        signif(min(object$eval.times), digits), "to",
        signif(max(object$eval.times), digits), "\n")
    cat("plot.times:", length(object$plot.times), "\n")
    if (length(object$plot.times) <= 10) {
        cat("plot.times values:", .show.times(object$plot.times), "\n")
    }
    cat("rv.times:", .show.times(object$rv.times), "\n")
    if (isTRUE(object$rmst)) {
        cat("fit.times.rmst:", .show.times(object$fit.times.rmst), "\n")
        cat("fit.times.rmst source:", object$time.info$fit.times.rmst.source, "\n")
        if (!is.null(object$time.info$rmst.message)) {
            cat("RMST note:", object$time.info$rmst.message, "\n")
        }
    }
    cat("fit.times source:", object$time.info$fit.times.source, "\n")
    cat("upper time source:", object$time.info$upper.time.source, "\n")
    cat("upper time:", signif(object$time.info$upper.time, digits), "\n")
    if (!is.null(object$time.info$G.cutoff)) {
        cat("G.cutoff:", object$time.info$G.cutoff, "\n")
    }
    invisible(object)
}

#' Plot default analysis time support
#'
#' @param x An object returned by \code{\link{npsa_times}()}.
#' @param ... Additional arguments.
#'
#' @return A \code{ggplot} object.
#'
#' @export
plot.npsa_times <- function(x, ...) {
    if (is.null(x$censor.df)) {
        stop("`x` does not contain censoring-support information.")
    }
    df <- x$censor.df
    p <- ggplot2::ggplot(df, ggplot2::aes(x = time, y = G.hat)) +
        ggplot2::geom_line(color = "black") +
        ggplot2::geom_hline(yintercept = x$time.info$G.cutoff,
                            linetype = "dashed", color = "red") +
        ggplot2::geom_vline(xintercept = max(x$fit.times),
                            linetype = "dotted", color = "blue") +
        ggplot2::xlab("Time") +
        ggplot2::ylab("Estimated censoring survival") +
        ggplot2::theme_bw()
    return(p)
}

.get.time.info <- function(time, event, fit.times, nuisance.options,
                           max.fit.times = 50, max.eval.times = 200,
                           G.cutoff = 0.05, verbose = FALSE) {
    if (is.null(nuisance.options)) nuisance.options <- list()
    time.info <- list()
    censor.df <- NULL

    if (is.null(fit.times)) {
        if (sum(event == 1) == 0) {
            stop("No uncensored events; cannot choose default fit.times.")
        }

        max.event <- max(time[event == 1])
        all.times <- sort(unique(time[time > 0 & time < max.event]))
        if (length(all.times) == 0) {
            all.times <- sort(unique(time[time > 0 & time <= max.event]))
        }
        if (length(all.times) == 0) {
            stop("No positive follow-up times available for default fit.times.")
        }
        all.times.full <- all.times

        if (sum(event == 0) == 0) {
            G.hat <- rep(1, length(all.times))
        } else {
            G.fit <- tryCatch({
                survival::survfit(survival::Surv(time, 1 - event) ~ 1)
            }, error = function(e) {
                NULL
            })
            G.hat <- tryCatch({
                summary(G.fit, times = all.times, extend = TRUE)$surv
            }, error = function(e) {
                rep(NA_real_, length(all.times))
            })
        }

        keep.idx <- which(is.finite(G.hat) & G.hat >= G.cutoff)
        upper.type <- "G.cutoff"
        if (length(keep.idx) == 0) {
            event.times <- sort(unique(time[event == 1 & time > 0]))
            upper.time <- as.numeric(stats::quantile(event.times, probs = 0.95,
                                                     type = 1, na.rm = TRUE))
            upper.type <- "event.quantile"
        } else {
            upper.time <- max(all.times[keep.idx])
        }
        if (!is.finite(upper.time)) upper.time <- max(all.times)

        censor.df <- data.frame(time = all.times.full,
                                G.hat = as.numeric(G.hat),
                                keep = all.times.full <= upper.time)
        all.times <- all.times[all.times <= upper.time]
        if (length(all.times) == 0) {
            all.times <- min(sort(unique(time[time > 0])))
        }
        fit.times <- all.times
        if (length(fit.times) > max.fit.times) {
            idx <- unique(round(seq(1, length(fit.times), length.out = max.fit.times)))
            fit.times <- fit.times[idx]
        }

        time.info$fit.times.source <- "automatic"
        time.info$upper.time.source <- upper.type
        time.info$upper.time <- max(fit.times)
        time.info$G.cutoff <- G.cutoff
        time.info$n.all.times <- length(all.times)
        if (verbose) {
            message("Using automatic fit.times with ", length(fit.times),
                    " time points up to ", signif(max(fit.times), 4), ".")
        }
    } else {
        if (!is.numeric(fit.times) || any(!is.finite(fit.times))) {
            stop("`fit.times` must be NULL or a finite numeric vector.")
        }
        fit.times <- sort(unique(fit.times))
        if (any(fit.times <= 0)) {
            fit.times <- fit.times[fit.times > 0]
            message("fit.times <= 0 removed.")
        }
        if (sum(event == 1) == 0) stop("No uncensored events; cannot perform estimation.")
        if (any(fit.times > max(time[event == 1]))) {
            fit.times <- fit.times[fit.times <= max(time[event == 1])]
            message("fit.times > max(time[event == 1]) removed.")
        }
        if (length(fit.times) == 0) {
            stop("No `fit.times` remain within the observed event-time range.")
        }
        time.info$fit.times.source <- "user"
        time.info$upper.time.source <- "user"
        time.info$upper.time <- max(fit.times)
        time.info$G.cutoff <- G.cutoff
        time.info$n.all.times <- length(fit.times)
    }

    if (is.null(censor.df) && sum(event == 1) > 0) {
        max.event <- max(time[event == 1])
        all.times <- sort(unique(time[time > 0 & time < max.event]))
        if (length(all.times) == 0) {
            all.times <- sort(unique(time[time > 0 & time <= max.event]))
        }
        if (length(all.times) > 0) {
            if (sum(event == 0) == 0) {
                G.hat <- rep(1, length(all.times))
            } else {
                G.fit <- tryCatch({
                    survival::survfit(survival::Surv(time, 1 - event) ~ 1)
                }, error = function(e) {
                    NULL
                })
                G.hat <- tryCatch({
                    summary(G.fit, times = all.times, extend = TRUE)$surv
                }, error = function(e) {
                    rep(NA_real_, length(all.times))
                })
            }
            censor.df <- data.frame(time = all.times,
                                    G.hat = as.numeric(G.hat),
                                    keep = all.times <= max(fit.times))
        }
    }

    if (is.null(nuisance.options$eval.times)) {
        eval.times <- sort(unique(time[time > 0 & time <= max(fit.times)]))
        if (length(eval.times) > max.eval.times) {
            idx <- unique(round(seq(1, length(eval.times), length.out = max.eval.times)))
            eval.times <- eval.times[idx]
        }
        nuisance.options$eval.times <- sort(unique(c(0, eval.times, max(fit.times))))
        time.info$eval.times.source <- "automatic"
    } else {
        time.info$eval.times.source <- "user"
    }

    time.info$fit.times <- fit.times
    time.info$eval.times <- nuisance.options$eval.times
    time.info$n.fit.times <- length(fit.times)
    time.info$n.eval.times <- length(nuisance.options$eval.times)
    time.info$max.fit.times <- max.fit.times
    time.info$max.eval.times <- max.eval.times

    list(fit.times = fit.times,
         nuisance.options = nuisance.options,
         time.info = time.info,
         censor.df = censor.df)
}

.get.report.times <- function(times, fit.times, default.all = TRUE,
                              label = "plot.times") {
    if (is.null(times)) {
        if (default.all) return(fit.times)
        return(numeric(0))
    }
    if (!is.numeric(times) || any(!is.finite(times))) {
        stop("`", label, "` must be NULL or a finite numeric vector.")
    }
    times <- sort(unique(times[times >= min(fit.times) & times <= max(fit.times)]))
    if (length(times) == 0) {
        stop("No `", label, "` remain within the fitted time range.")
    }
    matched <- sapply(times, function(t0) {
        fit.times[min(which(fit.times >= t0))]
    })
    matched <- sort(unique(as.numeric(matched)))
    if (!all(times %in% fit.times)) {
        message(label, " not in fit.times were matched to the nearest later fitted time.")
    }
    matched
}

#' Compute Effect Bounds (Internal)
#'
#' Internal utility function to compute the sensitivity analysis effect bounds.
#'
#' @param fit.times Numeric vector of times where bounds are evaluated.
#' @param theta.obs Numeric vector of observed differences.
#' @param psi Numeric vector of treatment assignment probabilities.
#' @param tau Truncation parameter.
#' @param sens.out Sensitivity parameter for U-T.
#' @param sens.trt Sensitivity parameter for U-A.
#' @param rho Correlation parameter.
#'
#' @return A list of lower and upper bounds at specified times.
#'
#' @keywords internal
.get.effect.bounds <- function(fit.times, theta.obs, psi, tau, sens.out, sens.trt, rho=1){
  effect.lower <- theta.obs - abs(rho)*sqrt(psi)*sqrt(tau)*sqrt(sens.out)*sqrt(sens.trt)
  effect.upper <- theta.obs + abs(rho)*sqrt(psi)*sqrt(tau)*sqrt(sens.out)*sqrt(sens.trt)

  res <- list(times=fit.times, effect.lower=effect.lower, effect.upper=effect.upper,
              sens.out=sens.out, sens.trt=sens.trt)
  return(res)
}

#' Compute Confidence Intervals for Effect Bounds (Internal)
#'
#' Internal utility function to compute pointwise and uniform confidence intervals for the estimated effect bounds.
#'
#' @param effect.bounds List containing lower and upper bounds for the effect.
#' @param psi Numeric vector of treatment probabilities at each time point.
#' @param tau Truncation parameter for survival times.
#' @param IF.vals.theta.obs Influence function values for observed treatment effect estimates.
#' @param IF.vals.psi Influence function values for treatment probability estimates.
#' @param IF.vals.tau Influence function values for tau.
#' @param conf.level Confidence level for intervals (default is 0.975).
#'
#' @return A list with confidence intervals for effect bounds, pointwise and uniform bands.
#'
#' @keywords internal
.bounds.confints <- function(effect.bounds, psi, tau,
                             IF.vals.theta.obs, IF.vals.psi, IF.vals.tau,
                             rho=1, band.end.pts = c(0,Inf), conf.level=.95, boot=10000, scale=TRUE){
  # check input
  # times, bias.lower, bias.upper,
  # IF.vals.theta.obs, IF.vals.psi
  # sens.out length should match with each other
  fit.times <- effect.bounds$times
  effect.lower <- effect.bounds$effect.lower
  effect.upper <- effect.bounds$effect.upper
  sens.out <- effect.bounds$sens.out
  sens.trt <- effect.bounds$sens.trt

  if (length(psi) != length(fit.times)) {
      stop("The length of `psi` must match the number of effect bound times.")
  }
  if (any(!is.finite(psi) | psi <= 0)) {
      stop("`psi` must contain finite positive values.")
  }
  if (length(tau) != 1 || !is.finite(tau) || tau <= 0) {
      stop("`tau` must be a finite positive value.")
  }

  n <- dim(IF.vals.theta.obs)[1]

  # IF for effect bound
  # IF.vals.effect.lower <- IF.vals.theta.obs # n*t
  # IF.vals.effect.upper <- IF.vals.theta.obs # n*t
  # inner.func.1 <- abs(rho)*sqrt(sens.out)*sqrt(sens.trt)*(1/2)/(sqrt(psi)*sqrt(tau)) # 1*t
  # for(j in 1:length(fit.times)) {
  #     inner.func.2 <- (tau*IF.vals.psi[,j]+psi[j]*IF.vals.tau)*inner.func.1[j]
  #     IF.vals.effect.lower[,j] <- IF.vals.theta.obs[,j] - inner.func.2
  #     IF.vals.effect.upper[,j] <- IF.vals.theta.obs[,j] + inner.func.2
  # }

  inner.func.1 <- abs(rho)*sqrt(sens.out)*sqrt(sens.trt)*(1/2)/(sqrt(psi)*sqrt(tau)) # 1*t
  inner.func.2 <- sweep((tau*IF.vals.psi+IF.vals.tau%*%rbind(psi)), 2, inner.func.1, '*') # n*t dot product 1*t column-wise multiplication
  IF.vals.effect.lower <- IF.vals.theta.obs - inner.func.2
  IF.vals.effect.upper <- IF.vals.theta.obs + inner.func.2

  # There might be NA in IF.gamma for RMST due to integration
  if (sum(is.na(IF.vals.effect.lower)) > 0) {
      row.notna.idx <- apply(IF.vals.effect.lower, 1, function(row) sum(is.na(row)) == 0)
      IF.vals.effect.lower <- IF.vals.effect.lower[row.notna.idx,]
      IF.vals.effect.upper <- IF.vals.effect.upper[row.notna.idx,]
      n <- nrow(IF.vals.effect.lower)
      message(sprintf("NA detected in IF.vals. Removed %d rows with missing values. Now N = %d.", sum(!row.notna.idx), n))
      if (n == 0) stop("No rows remain after removing missing IF values.")
  }

  # pointwise confidence intervals for effect bounds - correlated
  sigma2.l <- colMeans(IF.vals.effect.lower^2)
  sigma2.ul <- colMeans(IF.vals.effect.lower*IF.vals.effect.upper)
  sigma2.u <- colMeans(IF.vals.effect.upper^2)

  library(mvtnorm)
  c_alpha <- numeric(length(fit.times))
  for (t in 1:length(fit.times)) {
      cov.matrix <- matrix(c(sigma2.l[t], sigma2.ul[t], sigma2.ul[t], sigma2.u[t]),
                           nrow = 2, byrow = TRUE)
      epsilon <- rmvnorm(n=boot, mean=rep(0, 2), sigma = cov.matrix)
      epsilon[,2] <- - epsilon[,2]
      c_alpha[t] <- unname(quantile(apply(epsilon, 1, max), conf.level))
  }

  ptwise.bounds.lower <- effect.lower - c_alpha / sqrt(n)
  ptwise.bounds.upper <- effect.upper + c_alpha / sqrt(n)

  # pointwise transformation - correlated

  sigma2.trans.l <- (2/(1-effect.lower^2))^2*sigma2.l
  sigma2.trans.ul <- (4/((1-effect.lower^2)*(1-effect.upper^2)))*sigma2.ul
  sigma2.trans.u <- (2/(1-effect.upper^2))^2*sigma2.u

  trans.log <- function(x) log(1+x) - log(1-x)
  trans.log.inv <- function(x) (exp(x)-1)/(exp(x)+1)

  # scale to match uniform band transformation

  if (!scale){
      c_alpha.trans <- numeric(length(fit.times))
      for (t in 1:length(fit.times)) {
          cov.matrix <- matrix(c(sigma2.trans.l[t], sigma2.trans.ul[t], sigma2.trans.ul[t], sigma2.trans.u[t]),
                               nrow = 2, byrow = TRUE)
          epsilon <- mvtnorm::rmvnorm(n=boot, mean=rep(0, 2), sigma = cov.matrix)
          epsilon[,2] <- - epsilon[,2]
          c_alpha.trans[t] <- quantile(apply(epsilon, 1, max), conf.level)
      }

      ptwise.trans.l <- trans.log.inv(trans.log(pmin(pmax(effect.lower, -1), 1)) - c_alpha.trans / sqrt(n))
      ptwise.trans.u <- trans.log.inv(trans.log(pmin(pmax(effect.upper, -1), 1)) + c_alpha.trans / sqrt(n))
  } else{
      sigma.trans.l <- sqrt(sigma2.trans.l)
      sigma.trans.u <- sqrt(sigma2.trans.u)

      corr.trans <- sigma2.trans.ul / sqrt(sigma2.trans.l * sigma2.trans.u)
      corr.trans <- pmin(pmax(corr.trans, -1), 1)

      c_alpha.trans <- numeric(length(fit.times))

      for (t in 1:length(fit.times)) {
          cov.matrix <- matrix(
              c(1, corr.trans[t],
                corr.trans[t], 1),
              nrow = 2,
              byrow = TRUE
          )
          epsilon <- mvtnorm::rmvnorm(n = boot, mean = rep(0, 2), sigma = cov.matrix)
          epsilon[, 2] <- -epsilon[, 2]
          c_alpha.trans[t] <- unname(quantile(apply(epsilon, 1, max), conf.level))
      }

      ptwise.trans.l <- trans.log.inv(
          trans.log(pmin(pmax(effect.lower, -1), 1)) -
              c_alpha.trans * sigma.trans.l / sqrt(n)
      )

      ptwise.trans.u <- trans.log.inv(
          trans.log(pmin(pmax(effect.upper, -1), 1)) +
              c_alpha.trans * sigma.trans.u / sqrt(n)
      )
  }

  # uniform confidence band - correlated
  cut.index <- (fit.times >= band.end.pts[1] & fit.times <= band.end.pts[2])

  epsilon <- .estimate.limit.dist.bound(IF.vals.l=IF.vals.effect.lower,
                                        IF.vals.u=IF.vals.effect.upper,
                                        cut.index=cut.index)
  epsilon.l <- epsilon$epsilon.l
  epsilon.u <- epsilon$epsilon.u

  dist.null.sens <- apply(cbind(apply(epsilon.l,1,max), apply(-epsilon.u,1,max)),1,max)

  q_n <- unname(quantile(dist.null.sens, conf.level))

  uniform.bounds.lower <- effect.lower - q_n / sqrt(n)
  uniform.bounds.upper <- effect.upper + q_n / sqrt(n)
  uniform.bounds.lower[fit.times < band.end.pts[1] | fit.times > band.end.pts[2]] <- NA
  uniform.bounds.upper[fit.times < band.end.pts[1] | fit.times > band.end.pts[2]] <- NA

  # uniform transformation - correlated
  IF.trans.l <- IF.vals.effect.lower
  IF.trans.u <- IF.vals.effect.upper
  for(j in 1:length(effect.lower)) {
      IF.trans.l[,j] <- IF.vals.effect.lower[,j] * (2/(1-effect.lower[j]^2))
      IF.trans.u[,j] <- IF.vals.effect.upper[,j] * (2/(1-effect.upper[j]^2))
  }
  se.trans.l <- sqrt(colMeans(IF.trans.l^2))
  se.trans.u <- sqrt(colMeans(IF.trans.u^2))

  epsilon <- .estimate.limit.dist.bound(IF.vals.l=scale(IF.trans.l),
                                        IF.vals.u=scale(IF.trans.u),
                                        cut.index=cut.index)
  epsilon.l <- epsilon$epsilon.l
  epsilon.u <- epsilon$epsilon.u

  dist.null.sens <- apply(cbind(apply(epsilon.l,1,max), apply(-epsilon.u,1,max)),1,max)

  q_n <- unname(quantile(dist.null.sens, conf.level))

  uniform.trans.l <- trans.log.inv(trans.log(pmin(pmax(effect.lower, -1), 1)) - q_n*se.trans.l / sqrt(n))
  uniform.trans.l[fit.times < band.end.pts[1] | fit.times > band.end.pts[2]] <- NA
  uniform.trans.u <- trans.log.inv(trans.log(pmin(pmax(effect.upper, -1), 1)) + q_n*se.trans.u / sqrt(n))
  uniform.trans.u[fit.times < band.end.pts[1] | fit.times > band.end.pts[2]] <- NA

  res <- list(times=fit.times, effect.lower=effect.lower, effect.upper=effect.upper,
              ptwise.bounds.lower=ptwise.bounds.lower, ptwise.bounds.upper=ptwise.bounds.upper,
              ptwise.trans.lower=ptwise.trans.l, ptwise.trans.upper=ptwise.trans.u,
              uniform.bounds.lower=uniform.bounds.lower, uniform.bounds.upper=uniform.bounds.upper,
              uniform.trans.lower=uniform.trans.l, uniform.trans.upper=uniform.trans.u,
              sens.out=sens.out, sens.trt=sens.trt)
  return(res)
}

#' Convert Bounds and Confidence Intervals to Data Frame (Internal)
#'
#' Internal utility function to organize effect bounds and confidence intervals
#' into a tidy data frame for plotting or downstream analysis.
#'
#' @param bounds.conf.int List containing lower/upper confidence intervals and uniform bands.
#' @param theta.obs Numeric vector of observed treatment effect estimates.
#' @param d Number of dropped confounders (for labeling, optional).
#' @param transform Logical; whether to transform the bounds to survival differences.
#' @param time.zero Logical; whether to add a zero starting point at time = 0.
#'
#' @return A \code{data.frame} ready for plotting or reporting.
#'
#' @keywords internal
bounds2df <- function(bounds.conf.int, theta.obs, d=NULL, transform=TRUE, time.zero=TRUE){
    if (is.null(d)){
        d=0
    }

    bounds.df <- data.frame(times = bounds.conf.int$times, d = d)

  if (transform) {
    bounds.df$uniform.trans.lower <- bounds.conf.int$uniform.trans.lower
    bounds.df$ptwise.trans.lower <- bounds.conf.int$ptwise.trans.lower
    bounds.df$effect.lower <- bounds.conf.int$effect.lower
    bounds.df$theta.obs <- theta.obs
    bounds.df$effect.upper <- bounds.conf.int$effect.upper
    bounds.df$ptwise.trans.upper <- bounds.conf.int$ptwise.trans.upper
    bounds.df$uniform.trans.upper <- bounds.conf.int$uniform.trans.upper
  } else {
    bounds.df$uniform.bounds.lower <- bounds.conf.int$uniform.bounds.lower
    bounds.df$ptwise.bounds.lower <- bounds.conf.int$ptwise.bounds.lower
    bounds.df$effect.lower <- bounds.conf.int$effect.lower
    bounds.df$theta.obs <- theta.obs
    bounds.df$effect.upper <- bounds.conf.int$effect.upper
    bounds.df$ptwise.bounds.upper <- bounds.conf.int$ptwise.bounds.upper
    bounds.df$uniform.bounds.upper <- bounds.conf.int$uniform.bounds.upper
  }

  if (time.zero) {
    new_row <- as.data.frame(list(0, d, 0, 0, 0, 0, 0, 0, 0))
    names(new_row) <- names(bounds.df)
    if (transform) {
      new_row$uniform.trans.lower <- NA
      new_row$uniform.trans.upper <- NA
    } else {
      new_row$uniform.bounds.lower <- NA
      new_row$uniform.bounds.upper <- NA
    }

    bounds.df <- rbind(bounds.df, new_row)
  }
  return(bounds.df)
}


##############
# testing
##############
#' Compute Uniform Robustness Value (Internal)
#'
#' Internal utility function to compute the Uniform Robustness Value (RV)
#' for sensitivity analysis across multiple time points.
#'
#' @param theta.obs Numeric vector of observed treatment effect estimates.
#' @param rho Correlation parameter.
#' @param theta Hypothesized effect value.
#' @param verbose Logical; if TRUE, print uniform RV messages.
#'
#' @return Numeric value representing the uniform robustness value.
#'
#' @keywords internal
.get.uniform.RV <- function(theta.obs, psi, tau,
                            IF.vals.theta.obs, IF.vals.psi, IF.vals.tau,
                            rho=1, theta=0, conf.level=.95, seed=6741,
                            verbose=TRUE){
  set.seed(seed)

  # p-value under observed data
  # only proceed when p<0.05
  n <- nrow(IF.vals.theta.obs)
  # test.stat <- n^(1/2)*sum(abs(theta.obs))
  epsilon <- .estimate.limit.dist(IF.vals = IF.vals.theta.obs)
  # dist.null <- apply(epsilon, 1, function(x) {max(abs(x))})
  # dist.null <- rowSums(abs(epsilon)) # integration
  # dist.null <- replicate(1e4, sum(abs(rbind(rt(n, df = n - 1)/sqrt(n)) %*% IF.vals.theta.obs)))


  test.stat <- n^(1/2)*max(abs(theta.obs - theta))
  dist.null <- apply(epsilon, 1, function(x) {max(abs(x))})

  pvalue <- mean(dist.null > test.stat)
  if (verbose) cat("The p-value under no unobserved confounding is:", pvalue, "\n")

  if (pvalue < 1-conf.level) {

    if (verbose) message("Proceed to the test under unobserved confounding...")
    .get.pvalue.sens <- function(x){
      # IF function
      sens.all <- (x/sqrt(1-x))*abs(rho)
      inner.func.1 <- sens.all*(1/2)/(sqrt(psi)*sqrt(tau)) # 1*t
      inner.func.2 <- sweep((tau*IF.vals.psi+IF.vals.tau%*%rbind(psi)), 2, inner.func.1, '*') # n*t dot product 1*t column-wise multiplication

      IF.vals.effect.lower <- IF.vals.theta.obs - inner.func.2
      IF.vals.effect.upper <- IF.vals.theta.obs + inner.func.2

      # under H_0
      epsilon <- .estimate.limit.dist.bound(IF.vals.l=IF.vals.effect.lower,
                                            IF.vals.u=IF.vals.effect.upper)
      epsilon.l <- epsilon$epsilon.l
      epsilon.u <- epsilon$epsilon.u

      dist.null.sens <- apply(cbind(apply(epsilon.l,1,max), apply(-epsilon.u,1,max)),1,max)

      test.stat.sens <- max(n^(1/2)*max((theta.obs-sens.all*(sqrt(psi)*sqrt(tau))-theta)),
                            n^(1/2)*max(-(theta.obs+sens.all*(sqrt(psi)*sqrt(tau))-theta)))

      q_n <- unname(quantile(dist.null.sens, conf.level))

      # p-value
      pvalue.sens <- mean(dist.null.sens > test.stat.sens)

      return(pvalue.sens)
    }

    pvalue.root <- function(x){
      pvalue.sens <- .get.pvalue.sens(x)
      return(pvalue.sens - (1- conf.level))
    }

    uniform.RV <- tryCatch({
      uniroot(pvalue.root, c(0.0001,0.9999), tol = 0.0001)$root
    }, error = function(e) {
      message("An error occurred: ", e$message)
      NA
    })
  } else {
    if (verbose) message("The null hypothesis that the observed effect is zero cannot be rejected. Sensitivity analysis will not proceed.")
    uniform.RV <- NA
  }

  return(uniform.RV)
}

.estimate.limit.dist.bound <- function(IF.vals.l, IF.vals.u,
                                       cut.index=NULL, boot=10000){
  n <- nrow(IF.vals.l)
  t <- ncol(IF.vals.l)
  if (!is.null(cut.index)){
      IF.vals.l <- IF.vals.l[,cut.index, drop = FALSE]
      IF.vals.u <- IF.vals.u[,cut.index, drop = FALSE]
      t <- ncol(IF.vals.l)
  }

  cov.matrix.l <- matrix(0, nrow = t, ncol = t)
  for (i in 1:n) {
    cov.matrix.l <- cov.matrix.l + outer(IF.vals.l[i,], IF.vals.l[i,])
  }
  # if (trans) cov.matrix.l <- cov.matrix.l/(sigma2.trans.l)

  cov.matrix.u <- matrix(0, nrow = t, ncol = t)
  for (i in 1:n) {
    cov.matrix.u <- cov.matrix.u + outer(IF.vals.u[i,], IF.vals.u[i,])
  }
  # if (trans) cov.matrix.u <- cov.matrix.u/(sigma2.trans.u)

  cov.matrix.l.u <- matrix(0, nrow = t, ncol = t)
  for (i in 1:n) {
    cov.matrix.l.u <- cov.matrix.l.u + outer(IF.vals.l[i,], IF.vals.u[i,])
  }
  cov.matrix.u.l <- t(cov.matrix.l.u)
  # if (trans) {
  #     cov.matrix.l.u <- cov.matrix.l.u/(sqrt(sigma2.trans.l*sigma2.trans.u))
  #     cov.matrix.u.l <- cov.matrix.u.l/(sqrt(sigma2.trans.u*sigma2.trans.l))
  # }


  cov.matrix <- rbind(cbind(cov.matrix.l, cov.matrix.l.u), cbind(cov.matrix.u.l, cov.matrix.u))
  cov.matrix <- cov.matrix / n

  library(mvtnorm)
  epsilon <- rmvnorm(n=boot, mean=rep(0, 2*t), sigma = cov.matrix)
  epsilon.l <- epsilon[,1:t, drop = FALSE]
  epsilon.u <- epsilon[,(t+1):(2*t), drop = FALSE]

  return(list(epsilon.l=epsilon.l, epsilon.u=epsilon.u))
}

.estimate.limit.dist <- function(IF.vals, boot=10000){

    n <- nrow(IF.vals)
    t <- ncol(IF.vals)

    cov.matrix <- matrix(0, nrow = t, ncol = t)
    for (i in 1:n) {
        cov.matrix <- cov.matrix + outer(IF.vals[i,], IF.vals[i,])
    }
    cov.matrix <- cov.matrix / n

    library(mvtnorm)
    epsilon <- rmvnorm(n=boot, mean=rep(0, t), sigma = cov.matrix)
    return(epsilon)
}

# for (x in seq(0.02,0.03,by=0.001)){
#     cat(x, .get.pvalue.sens(x), "\n")
# }

# validate
# inner.func.1 <- (x/sqrt(1-x))*abs(rho)*(1/2)/(sqrt(psi[3])*sqrt(tau))
# inner.func.2 <- (tau*IF.vals.psi[,3]+psi[3]*IF.vals.tau)*inner.func.1 # n*1
# (IF.vals.theta.obs[,3] - inner.func.2)[1:10]

#################
#' Compute Pointwise Robustness Value (Internal)
#'
#' Internal utility function to compute the Pointwise Robustness Value (RV)
#' for sensitivity analysis at specific evaluation times.
#'
#' @param t0 Numeric evaluation time.
#' @param fit.times Numeric vector of fitted times.
#' @param theta.obs Numeric vector of observed treatment effect estimates.
#' @param psi Numeric vector of estimated psi values.
#' @param tau Numeric estimated tau value.
#' @param IF.vals.theta.obs Influence function values for the observed effect.
#' @param IF.vals.psi Influence function values for psi.
#' @param IF.vals.tau Influence function values for tau.
#' @param rho Correlation parameter.
#' @param theta Hypothesized effect value.
#' @param conf.bounds Logical; if TRUE, compute MIRV using confidence bounds.
#' @param transform Logical; whether to use transformed pointwise MIRV calculation.
#' @param conf.level Confidence level.
#' @param verbose Logical; if TRUE, print the no-sensitivity confidence interval.
#' @param boot Number of simulations for the pointwise confidence calculation.
#'
#' @return A numeric vector of pointwise robustness values corresponding to each evaluation time.
#'
#' @keywords internal
.get.RV <- function(t0, fit.times, theta.obs, psi, tau,
                    IF.vals.theta.obs, IF.vals.psi, IF.vals.tau,
                    rho=1, theta=0, conf.bounds=TRUE, transform=FALSE, conf.level=.95, verbose=TRUE, boot=10000){
  res.RV <- NULL

  k <- min(which(fit.times >= t0))
  theta.obs.t0 = theta.obs[k]
  psi.t0 = psi[k]

  v <- (theta.obs.t0 - theta)^2/(rho^2*psi.t0*tau)
  bounds.RV <- (- v + sqrt(v^2 + 4*v))/2

  res.RV <- list(t0=t0, bounds.RV=bounds.RV, theta=theta, bounds.int.RV=NA, conf.level=conf.level)

  if (conf.bounds){
    n <- dim(IF.vals.theta.obs)[1]
    epsilon.base <- matrix(stats::rnorm(boot * n), nrow = boot, ncol = n)

    # pointwise confidence intervals as function of sensitivity parameters
    bounds.senspar <- function(x, lower.b=TRUE){
      inner.func.1 <- (x/sqrt(1-x))*abs(rho)*(1/2)/(sqrt(psi.t0)*sqrt(tau))
      inner.func.2 <- (tau*IF.vals.psi[,k]+psi.t0*IF.vals.tau)*inner.func.1 # n*1

      IF.vals.effect.lower <- IF.vals.theta.obs[,k] - inner.func.2
      IF.vals.effect.upper <- IF.vals.theta.obs[,k] + inner.func.2

      effect.lower.sp <- theta.obs.t0 - abs(rho)*sqrt(psi.t0)*sqrt(tau)*(x/sqrt(1-x))
      effect.upper.sp <- theta.obs.t0 + abs(rho)*sqrt(psi.t0)*sqrt(tau)*(x/sqrt(1-x))

      if (!transform){
          IF.vals.effect <- cbind(IF.vals.effect.lower, IF.vals.effect.upper)
      } else {
          IF.vals.effect <- cbind(IF.vals.effect.lower * (2/(1-effect.lower.sp^2)),
                                  IF.vals.effect.upper * (2/(1-effect.upper.sp^2)))

          trans.log <- function(x) log(1+x) - log(1-x)
          trans.log.inv <- function(x) (exp(x)-1)/(exp(x)+1)
      }

      epsilon <- (epsilon.base / sqrt(n)) %*% IF.vals.effect
      epsilon[,2] <- - epsilon[,2]
      c_alpha <- unname(quantile(apply(epsilon, 1, max), conf.level))

      if (lower.b) {
          if (transform) {
              val <- trans.log(pmin(pmax(effect.lower.sp, -1), 1)) - c_alpha / sqrt(n)
              return(trans.log.inv(val) - theta)
          } else {
              return(effect.lower.sp - c_alpha / sqrt(n) - theta)
          }
      } else {
          if (transform) {
              val <- trans.log(pmin(pmax(effect.upper.sp, -1), 1)) + c_alpha / sqrt(n)
              return(trans.log.inv(val) - theta)
          } else {
              return(effect.upper.sp + c_alpha / sqrt(n) - theta)
          }
      }
    }

    # confidence intervals for theta.obs.t0 when sensitivity parameters are 0
    l.sp.0 <- bounds.senspar(0, lower.b = TRUE)
    u.sp.0 <- bounds.senspar(0, lower.b = FALSE)
    if (verbose) cat("pointwise confidence interval for theta.obs.t0 when sensitivity parameters are 0:", c(l.sp.0, u.sp.0), "\n")

    if (l.sp.0 <= 0 & u.sp.0 >= 0) {
      if (verbose) message("Pointwise CI cover the hypothesized value of theta; robustness values calculation for the lower/upper limit unnecessary.")
      bounds.int.RV <- 0
      lower.b <- NULL
    } else {
      lower.b <- ifelse(l.sp.0 > 0, TRUE, FALSE)
      bounds.int.RV <- tryCatch({
        x.grid <- unique(c(0, seq(0.0001, 0.99, length.out = 200)))
        bounds.vals <- sapply(x.grid, function(x) {
          tryCatch({
            out <- bounds.senspar(x, lower.b = lower.b)
            ifelse(length(out) == 1 && is.finite(out), out, NA)
          }, error = function(e) {
            NA
          })
        })
        bounds.vals[!is.finite(bounds.vals)] <- NA

        change.idx <- which(
          !is.na(bounds.vals[-length(bounds.vals)]) &
            !is.na(bounds.vals[-1]) &
            bounds.vals[-length(bounds.vals)] * bounds.vals[-1] <= 0
        )

        if (length(change.idx) == 0) {
          if (verbose) message("MIRV could not be computed: no finite sensitivity value in [0, 0.99] moved the pointwise confidence bound to the hypothesized effect.")
          NA
        } else {
          idx <- change.idx[1]
          if (bounds.vals[idx] == 0) {
            x.grid[idx]
          } else if (bounds.vals[idx + 1] == 0) {
            x.grid[idx + 1]
          } else {
            uniroot(bounds.senspar, c(x.grid[idx], x.grid[idx + 1]),
                    tol = 0.0001, lower.b = lower.b,
                    f.lower = bounds.vals[idx],
                    f.upper = bounds.vals[idx + 1])$root
          }
        }
      }, error = function(e) {
        if (verbose) message("An error occurred: ", e$message)
        NA
      })
    }

    res.RV <- list(t0=t0, bounds.RV=bounds.RV, theta=theta, bounds.int.RV=bounds.int.RV, conf.level=conf.level, lower.b=lower.b)
  }

  return(res.RV)
}
