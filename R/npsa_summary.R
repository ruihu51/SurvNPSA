#' Report Robustness Values
#'
#' @param transform Logical; whether to use transformed pointwise MIRV calculation.
#'
#' @keywords internal
.report.RV <- function(rv.times, result, rho = 1, theta = 0,
                       conf.level = .95, transform = FALSE,
                       unif = TRUE, t.lower, t.upper) {
  res.list <- list()

  if (any(rv.times > max(result$fit.times))) {
      message("Some rv.times > maximum observed event time - removed for RV computation.")
      rv.times <- rv.times[rv.times <= max(result$fit.times)]
  }

  for (t0 in rv.times){
    res.RV <- .get.RV(
      t0,
      fit.times = result$fit.times,
      theta.obs = result$obs.comps.df$theta.obs,
      psi = result$obs.comps.df$psi,
      tau = result$tau,
      IF.vals.theta.obs = result$IF.vals.theta.obs,
      IF.vals.psi = result$IF.vals.psi,
      IF.vals.tau = result$IF.vals.tau,
      rho = rho,
      theta = theta,
      transform = transform,
      conf.level = conf.level
    )

    res.list[[length(res.list) + 1]] <- list(
      t0 = res.RV$t0,
      theta = res.RV$theta,
      RV = res.RV$bounds.RV,
      MIRV = res.RV$bounds.int.RV,
      conf.level = res.RV$conf.level,
      lower.b = if (is.null(res.RV$lower.b)) NA else res.RV$lower.b,
      rho = rho
    )
  }

  res.table <- do.call(rbind, lapply(res.list, as.data.frame))
  rownames(res.table) <- NULL

  out <- list(res.table = res.table)

  if (unif){
    unif.idx <- which(result$fit.times >= t.lower & result$fit.times <= t.upper)
    unif.RV <- .get.uniform.RV(
      theta.obs = result$obs.comps.df$theta.obs[unif.idx],
      psi = result$obs.comps.df$psi[unif.idx],
      tau = result$tau,
      IF.vals.theta.obs = result$IF.vals.theta.obs[,unif.idx],
      IF.vals.psi = result$IF.vals.psi[,unif.idx],
      IF.vals.tau = result$IF.vals.tau,
      rho = rho,
      theta = theta,
      conf.level = conf.level
    )
    out$unif.RV <- unif.RV
    out$unif.idx <- unif.idx
  }

  class(out) <- "reportRV"
  return(out)
}

#' Summarize Reported Robustness Value (RV)
#'
#' Summary method for objects of class \code{reportRV}.
#'
#' @param object An object of class \code{reportRV}.
#' @param digits Number of digits for printing.
#' @param ... Additional arguments (currently unused).
#'
#' @return Printed robustness value results.
#'
#' @export
#' @method summary reportRV
summary.reportRV <- function(object, digits = 3, ...) {
  cat("Robustness Value Report\n")
  cat("------------------------\n")

  tbl <- object$res.table
  num.cols <- sapply(tbl, is.numeric)

  tbl[, num.cols] <- lapply(tbl[, num.cols, drop = FALSE], function(x) round(x, digits))

  print(tbl, row.names = FALSE)

  cat("\nFootnote:\n")
  cat("\u00B9 MIRV = 0 indicates that the pointwise confidence interval already covers the hypothesized value of theta; robustness value calculation for the lower/upper limit is unnecessary.\n")

  if (!is.null(object$unif.RV)) {
    cat("\nUniform Robustness Value available.\n")
  }
}

#' Interpret Robustness Values
#'
#' @param object An object returned by \code{npsa_surv()}.
#' @param t0 Time point for pointwise RV or MIRV interpretation.
#' @param type Which robustness value to interpret. Use \code{"RV"},
#'   \code{"MIRV"}, or \code{"URV"}.
#' @param var_names Optional confounder names. If \code{NULL}, names stored in
#'   \code{object} are used.
#'
#' @return An object of class \code{interpretRV}.
#'
#' @export
interpret.RV <- function(object, t0 = NULL, type = c("RV", "MIRV", "URV"),
                         var_names = NULL) {
    type <- match.arg(type)

    if (!inherits(object, "npsa_surv")) {
        stop("`object` must be an object returned by `npsa_surv()`.")
    }
    if (is.null(object$res.RV)) {
        stop("`object` does not contain RV results. Please run `npsa_surv()` with `rv.options$rv.times`.")
    }
    if (is.null(object$senspar.df$sens.df) || is.null(object$senspar.df$sens.df.mean)) {
        stop("`object` does not contain sensitivity parameter results for RV interpretation.")
    }

    if (is.null(var_names)) var_names <- object$var_names
    if (is.null(var_names)) {
        stop("Please provide `var_names`, or run `npsa_surv()` with confounder names.")
    }

    if (type == "URV") {
        return(.interpret.URV(object$result$fit.times, object$res.RV,
                              object$senspar.df$sens.df,
                              object$senspar.df$sens.df.mean,
                              var_names))
    }

    if (is.null(t0)) stop("Please provide `t0` for RV or MIRV interpretation.")

    .interpret.RV(t0 = t0, res.RV = object$res.RV,
                  sens.df = object$senspar.df$sens.df,
                  sens.df.mean = object$senspar.df$sens.df.mean,
                  var_names = var_names, type = type)
}

#' Interpret Pointwise Robustness Values
#'
#' @keywords internal
.interpret.RV <- function(t0, res.RV, sens.df, sens.df.mean, var_names,
                          type = c("RV", "MIRV")) {
    type <- match.arg(type)

    res.table <- res.RV$res.table

    if (!(t0 %in% res.table$t0)) {
        stop("No robustness values result for time t0.")
    }

    n_var <- length(var_names)
    rv <- res.table[res.table$t0 == t0, type]
    rv <- as.numeric(rv[1])

    if (is.na(rv)) {
        stop(sprintf("No %s result is available at time t0.", type))
    }

    sp.point <- rv^2 / (1 - rv)
    half_d <- ceiling(n_var * 0.5)

    sens.t0 <- sens.df %>%
        mutate(sens.par = C.Y.sq * C.A.sq,
               confounder = var_names[j]) %>%
        filter(near(t, t0))

    sens.mean.t0 <- sens.df.mean %>%
        filter(near(t, t0)) %>%
        arrange(d) %>%
        mutate(sig.point = sens.par > sp.point)

    out.1 <- sens.t0 %>%
        filter(d == 1, sens.par > sp.point) %>%
        pull(confounder)
    out.1 <- if (length(out.1) == 0) NULL else out.1

    change.idx <- which(diff(sens.mean.t0$sig.point) == 1)
    out.d <- if (length(change.idx) > 0) {
        c(sens.mean.t0$d[change.idx], sens.mean.t0$d[change.idx + 1])
    } else {
        NULL
    }

    out <- mean(
        sens.t0 %>%
            filter(d == half_d) %>%
            mutate(value = sens.par <= sp.point) %>%
            pull(value)
    )

    out.half <- if (is.na(out) || out == 1) {
        NULL
    } else {
        out
    }

    # ------ FINAL summary output --------

    summary_table <- tibble::tibble(
        Method = c("Leave-one-out", "Leave-d-out", "Leave-half-out"),
        Interpretation = c(
            if (is.null(out.1)) "None" else paste(out.1, collapse = ", "),
            if (is.null(out.d)) "None" else paste0("d=", paste(out.d, collapse = " and d=")),
            if (is.null(out.half)) "None" else paste0(round((out.half) * 100, 1), "th percentile")
        )
    )

    title_line <- paste0("Interpretation of ", type, " at time $t=", t0, "$")

    out <- list(
        title = title_line,
        table = summary_table,
        type = type,
        t0 = t0,
        rv = rv,
        sens.par = sp.point
    )
    class(out) <- "interpretRV"
    return(out)
}


.interpret.URV <- function(fit.times, res.RV, sens.df, sens.df.mean, var_names){

    unif.RV <- res.RV$unif.RV
    unif.idx <- res.RV$unif.idx
    sp.unif <- unif.RV^2 / (1 - unif.RV)

    n_var <- length(var_names)
    out.1 <- out.d <- out.half <- NULL

    out.unif <- sens.df %>%
        filter(
            t >= min(fit.times[unif.idx], na.rm = TRUE) &
            t <= max(fit.times[unif.idx], na.rm = TRUE)
        ) %>%
        mutate(sens.par = C.Y.sq * C.A.sq) %>%
        group_by(d, j) %>%
        summarize(sens.par = max(sens.par), .groups = "drop") %>%
        ungroup() %>%
        filter(d == 1, sens.par > sp.unif) %>%
        mutate(confounder = var_names[j])

    out.1 <- if (nrow(out.unif) == 0) NULL else out.unif$confounder

    out.unif <- sens.df %>%
        filter(
            t >= min(fit.times[unif.idx], na.rm = TRUE) &
            t <= max(fit.times[unif.idx], na.rm = TRUE)
        ) %>%
        mutate(sens.par = C.Y.sq * C.A.sq) %>%
        group_by(d, j) %>%
        summarize(sens.par = max(sens.par), .groups = "drop") %>%
        ungroup() %>%
        group_by(d) %>%
        summarize(sens.par = mean(sens.par), .groups = "drop") %>%
        arrange(d) %>%
        mutate(sig.unif = sens.par > sp.unif)

    change.idx <- which(diff(out.unif$sig.unif) == 1)
    out.d <- if (length(change.idx) > 0) c(out.unif$d[change.idx], out.unif$d[change.idx + 1]) else NULL

    out <- mean(
        sens.df %>%
            filter(
                t >= min(fit.times[unif.idx], na.rm = TRUE) &
                t <= max(fit.times[unif.idx], na.rm = TRUE)
            ) %>%
            mutate(sens.par = C.Y.sq * C.A.sq) %>%
            group_by(d, j) %>%
            summarize(sens.par = max(sens.par), .groups = "drop") %>%
            ungroup() %>%
            filter(d == ceiling(n_var * 0.5)) %>%
            mutate(value = sens.par <= sp.unif) %>%
            pull(value)
    )

    out.half <- if (!is.null(out) && !is.na(out) && out == 1) {
        NULL
    } else {
        out
    }

    # ------ FINAL summary output --------

    summary_table <- tibble::tibble(
        Method = c("Leave-one-out", "Leave-d-out", "Leave-half-out"),
        Interpretation = c(
            if (is.null(out.1)) "None" else paste(out.1, collapse = ", "),
            if (is.null(out.d)) "None" else paste0("d=", paste(out.d, collapse = " and d=")),
            if (is.null(out.half)) "None" else paste0(round((out.half) * 100, 1), "th percentile")
        )
    )

    title_line <- paste0("Interpretation of URV")

    out <- list(
        title = title_line,
        table = summary_table
    )
    class(out) <- "interpretRV"
    return(out)

}

#' Summarize Robustness Value (RV) Results
#'
#' Summary method for objects of class \code{interpretRV}, typically created
#' by sensitivity analysis functions in the \code{npsaSurv} package.
#'
#' @param object An object of class \code{interpretRV}.
#' @param ... Additional arguments (currently unused).
#'
#' @return A printed summary of robustness values and key interpretations.
#'
#' @export
#' @method summary interpretRV
summary.interpretRV <- function(object, ...) {
    cat(object$title, "\n\n")
    print(object$table)
    invisible(object)
}

#' Estimate observed components
#'
#' @keywords internal
.report.bounds <- function(plot.times, result, rho=1, band.end.pts = c(0,Inf), conf.level=.95, boot=10000,
                           sens.df.mean = NULL, num_drop = NULL, pct_drop = NULL, n_var = NULL,
                           rmst = TRUE, sens.rmst.df.mean = NULL, transform = TRUE, scale = TRUE) {

    if (any(plot.times > max(result$fit.times))) {
        message("Some plot.times > maximum observed event time - removed for plot.")
        plot.times <- plot.times[plot.times <= max(result$fit.times)]
    }

    if (is.null(sens.df.mean)) {

        obs.est.idx <- sapply(plot.times, function(x) {
            which(near(x, result$fit.times))
        })

        # Estimate lower and upper bounds without sensitivity (assume zero sensitivity)
        effect.bounds <- .get.effect.bounds(
            fit.times = result$fit.times[obs.est.idx],
            theta.obs = result$obs.comps.df$theta.obs[obs.est.idx],
            psi = result$obs.comps.df$psi[obs.est.idx],
            tau = result$tau,
            sens.out = rep(0, length(plot.times)),
            sens.trt = 0,
            rho = rho
        )

        bounds.conf.int <- .bounds.confints(
            effect.bounds,
            psi = result$obs.comps.df$psi[obs.est.idx],
            tau = result$tau,
            IF.vals.theta.obs = result$IF.vals.theta.obs[, obs.est.idx],
            IF.vals.psi = result$IF.vals.psi[, obs.est.idx],
            IF.vals.tau = result$IF.vals.tau,
            rho = rho,
            band.end.pts = band.end.pts,
            conf.level = conf.level,
            scale = scale,
            boot = boot
        )

        bounds.df <- bounds2df(bounds.conf.int, theta.obs = result$obs.comps.df$theta.obs[obs.est.idx],
                               d = NULL, transform = transform, time.zero = TRUE)

        bounds.df.rmst <- NULL

        if (rmst) {

            effect.bounds.rmst <- .get.effect.bounds(
                fit.times = result$fit.times.rmst,
                theta.obs = result$rmst.obs,
                psi = result$gamma.est,
                tau = result$tau,
                sens.out = rep(0, length(result$rmst.obs)),
                sens.trt = 0,
                rho = rho
            )

            bounds.conf.int.rmst <- .bounds.confints(
                effect.bounds.rmst,
                psi = result$gamma.est,
                tau = result$tau,
                IF.vals.theta.obs = result$IF.vals.rmst.obs,
                IF.vals.psi = result$IF.vals.gamma,
                IF.vals.tau = result$IF.vals.tau,
                rho = rho,
                band.end.pts = band.end.pts,
                conf.level = conf.level,
                scale = scale,
                boot = boot
            )

            bounds.df.rmst <- bounds2df(bounds.conf.int.rmst, theta.obs = result$rmst.obs,
                                        d = NULL, transform = FALSE, time.zero = FALSE)
        }

        class(bounds.df) <- c("boundsdf", "data.frame")

        return(list(bounds.df = bounds.df, bounds.df.rmst = bounds.df.rmst))

    } else {

        if (is.null(num_drop) && is.null(pct_drop)) {
            stop("You must specify either 'num_drop' or 'pct_drop'.")
        }
        if (!is.null(num_drop) && !is.null(pct_drop)) {
            stop("Specify only one of 'num_drop' or 'pct_drop'.")
        }

        if (!is.null(pct_drop)) {
            num_drop <- unique(ceiling(pct_drop * n_var))
            num_drop <- num_drop[num_drop >= 1 & num_drop < n_var]
        }

        invalid_values <- setdiff(num_drop, sens.df.mean$d)
        if (length(invalid_values) > 0) {
            stop("Invalid `num_drop` value(s): ", paste(invalid_values, collapse = ", "), ". Not provided in senspar.")
        }

        bounds.df <- data.frame()
        bounds.df.rmst <- data.frame()

        for (d in num_drop) {

            sens.out.true.input <- as.vector(sens.df.mean[sens.df.mean$d == d, "sens.par"])$sens.par
            sens.trt.true <- 1

            senspar.idx <- sapply(plot.times, function(x) {
                which(near(x, sens.df.mean$t[sens.df.mean$d == d]))
            })

            obs.est.idx <- sapply(plot.times, function(x) {
                which(near(x, result$fit.times))
            })

            effect.bounds <- .get.effect.bounds(
                fit.times = result$fit.times[obs.est.idx],
                theta.obs = result$obs.comps.df$theta.obs[obs.est.idx],
                psi = result$obs.comps.df$psi[obs.est.idx],
                tau = result$tau,
                sens.out = sens.out.true.input[senspar.idx],
                sens.trt = sens.trt.true,
                rho = rho
            )

            bounds.conf.int <- .bounds.confints(
                effect.bounds,
                psi = result$obs.comps.df$psi[obs.est.idx],
                tau = result$tau,
                IF.vals.theta.obs = result$IF.vals.theta.obs[, obs.est.idx],
                IF.vals.psi = result$IF.vals.psi[, obs.est.idx],
                IF.vals.tau = result$IF.vals.tau,
                conf.level = conf.level,
                scale = scale
            )

            df <- bounds2df(bounds.conf.int, theta.obs = result$obs.comps.df$theta.obs[obs.est.idx],
                            d = d, transform = transform, time.zero = TRUE)

            bounds.df <- rbind(bounds.df, df)

            if (rmst) {
                if (is.null(sens.rmst.df.mean)) {
                    stop("You must provide `sens.rmst.df.mean` when `rmst = TRUE`.")
                }

                sens.out.true.input.rmst <- as.vector(sens.rmst.df.mean[sens.rmst.df.mean$d == d, "sens.par"])$sens.par

                effect.bounds.rmst <- .get.effect.bounds(
                    fit.times = result$fit.times.rmst,
                    theta.obs = result$rmst.obs,
                    psi = result$gamma.est,
                    tau = result$tau,
                    sens.out = sens.out.true.input.rmst,
                    sens.trt = 1,
                    rho = rho
                )

                bounds.conf.int.rmst <- .bounds.confints(
                    effect.bounds.rmst,
                    psi = result$gamma.est,
                    tau = result$tau,
                    IF.vals.theta.obs = result$IF.vals.rmst.obs,
                    IF.vals.psi = result$IF.vals.gamma,
                    IF.vals.tau = result$IF.vals.tau,
                    conf.level = conf.level,
                    scale = scale,
                )

                df.rmst <- bounds2df(bounds.conf.int.rmst, theta.obs = result$rmst.obs,
                                     d = d, transform = FALSE, time.zero = FALSE)

                bounds.df.rmst <- rbind(bounds.df.rmst, df.rmst)
            }

        }

        class(bounds.df) <- c("boundsdf", "data.frame")

        return(list(bounds.df = bounds.df, bounds.df.rmst = bounds.df.rmst))

    }
}

#' Plot boundsdf object
#'
#' Plot method for objects of class \code{boundsdf}.
#'
#' @param x An object of class \code{boundsdf}.
#' @param ... Additional arguments (currently unused).
#'
#' @return A \code{ggplot} object showing sensitivity bounds.
#'
#' @export
#' @method plot boundsdf
plot.boundsdf <- function(x, ...) {
    if ("ptwise.trans.lower" %in% names(x)) {
        x$ptwise.lower <- x$ptwise.trans.lower
        x$ptwise.upper <- x$ptwise.trans.upper
        x$uniform.lower <- x$uniform.trans.lower
        x$uniform.upper <- x$uniform.trans.upper
    } else {
        x$ptwise.lower <- x$ptwise.bounds.lower
        x$ptwise.upper <- x$ptwise.bounds.upper
        x$uniform.lower <- x$uniform.bounds.lower
        x$uniform.upper <- x$uniform.bounds.upper
    }

    x %>%
        mutate(setting = paste0("Drop ", d, " confounder", ifelse(d > 1, "s", ""))) %>%
        ggplot(aes(x = times)) +
        geom_line(aes(y = theta.obs, linetype = "Observed Effect", color = "Observed Effect")) +
        geom_line(aes(y = effect.lower, linetype = "Effect Bounds", color = "Effect Bounds")) +
        geom_line(aes(y = effect.upper, linetype = "Effect Bounds", color = "Effect Bounds")) +
        geom_line(aes(y = ptwise.lower, linetype = "Pointwise CI", color = "Pointwise CI")) +
        geom_line(aes(y = ptwise.upper, linetype = "Pointwise CI", color = "Pointwise CI")) +
        geom_line(aes(y = uniform.lower, linetype = "Uniform Bands", color = "Uniform Bands")) +
        geom_line(aes(y = uniform.upper, linetype = "Uniform Bands", color = "Uniform Bands")) +
        scale_color_manual(values = c(
            "Observed Effect" = "black",
            "Effect Bounds" = "red",
            "Pointwise CI" = "blue",
            "Uniform Bands" = "brown"
        )) +
        scale_linetype_manual(values = c(
            "Observed Effect" = "solid",
            "Effect Bounds" = "dashed",
            "Pointwise CI" = "dotdash",
            "Uniform Bands" = "longdash"
        )) +
        labs(linetype = "Type", color = "Type") +
        xlab("Time") +
        ylab("Survival difference (treatment - control)") +
        theme_bw() +
        theme(
            legend.position = "bottom",
            text = element_text(size = 12),
            legend.text = element_text(size = 12),
            legend.key.width = unit(0.8, "cm"),
            legend.title = element_blank(),
            panel.grid.minor = element_blank()
        ) +
        facet_wrap(~setting, scales = "free_y")
}
