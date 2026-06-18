#' Report Robustness Values
#'
#' @param transform Logical; whether to use transformed pointwise MIRV calculation.
#' @param verbose Logical; if TRUE, print pointwise and uniform RV messages.
#'
#' @keywords internal
.report.RV <- function(rv.times, result, rho = 1, theta = 0,
                       conf.level = .95, transform = FALSE,
                       verbose = FALSE, unif = TRUE, t.lower, t.upper) {
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
      conf.level = conf.level,
      verbose = verbose
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
      conf.level = conf.level,
      verbose = verbose
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
        stop("`object` does not contain RV results. Please run `npsa_surv()` with default RV times or provide `rv.options$rv.times`.")
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

.npsa_sp_from_rv <- function(x) {
    out <- x^2 / (1 - x)
    out[is.nan(out)] <- NA
    out
}

.npsa_summary_time_rows <- function(df, report.times, time.col = "times") {
    if (is.null(df) || nrow(df) == 0) return(df)
    if (is.null(report.times) || length(report.times) == 0) return(df)

    out <- data.frame()
    if ("d" %in% names(df)) {
        split.df <- split(df, df$d)
    } else {
        split.df <- list(df)
    }

    for (one.df in split.df) {
        if (nrow(one.df) == 0) next
        idx <- sapply(report.times, function(t0) {
            which.min(abs(one.df[[time.col]] - t0))
        })
        idx <- unique(idx)
        out <- rbind(out, one.df[idx, , drop = FALSE])
    }
    rownames(out) <- NULL
    out
}

.npsa_bounds_summary_table <- function(df, report.times = NULL,
                                       effect.name = "surv.diff",
                                       include.d = TRUE) {
    if (is.null(df) || nrow(df) == 0) return(NULL)
    df <- as.data.frame(df)
    df <- .npsa_summary_time_rows(df, report.times, time.col = "times")
    if (is.null(df) || nrow(df) == 0) return(NULL)

    if ("ptwise.trans.lower" %in% names(df)) {
        ptwise.lower <- df$ptwise.trans.lower
        ptwise.upper <- df$ptwise.trans.upper
        unif.lower <- df$uniform.trans.lower
        unif.upper <- df$uniform.trans.upper
    } else {
        ptwise.lower <- df$ptwise.bounds.lower
        ptwise.upper <- df$ptwise.bounds.upper
        unif.lower <- df$uniform.bounds.lower
        unif.upper <- df$uniform.bounds.upper
    }

    out <- data.frame(time = df$times)
    if (include.d && "d" %in% names(df)) out$d <- df$d
    out[[effect.name]] <- df$theta.obs
    out$lower.bound <- df$effect.lower
    out$upper.bound <- df$effect.upper
    out$ptwise.lower <- ptwise.lower
    out$ptwise.upper <- ptwise.upper
    out$unif.lower <- unif.lower
    out$unif.upper <- unif.upper
    out$ci.includes.0 <- out$ptwise.lower <= 0 & out$ptwise.upper >= 0

    if (effect.name == "rmst.diff") names(out)[names(out) == "time"] <- "rmst.time"
    rownames(out) <- NULL
    out
}

.npsa_senspar_summary_table <- function(senspar.df, report.times = NULL) {
    if (is.null(senspar.df) || is.null(senspar.df$sens.df.mean)) return(NULL)
    df <- as.data.frame(senspar.df$sens.df.mean)
    if (nrow(df) == 0) return(NULL)
    df <- .npsa_summary_time_rows(df, report.times, time.col = "t")

    drop.info <- NULL
    if (!is.null(senspar.df$drop.sets)) {
        drop.info <- as.data.frame(senspar.df$drop.sets)
    } else if (!is.null(senspar.df$sens.df) && "drop.name" %in% names(senspar.df$sens.df)) {
        drop.info <- unique(as.data.frame(senspar.df$sens.df[, c("j", "d", "drop.name")]))
    }

    n.drop.sets <- rep(NA_integer_, nrow(df))
    drop.examples <- rep(NA_character_, nrow(df))

    if (!is.null(drop.info) && nrow(drop.info) > 0) {
        for (i in seq_len(nrow(df))) {
            tmp <- drop.info[drop.info$d == df$d[i], , drop = FALSE]
            n.drop.sets[i] <- length(unique(tmp$j))
            if ("drop.name" %in% names(tmp)) {
                examples <- unique(tmp$drop.name)
            } else if ("drop.index" %in% names(tmp)) {
                examples <- unique(tmp$drop.index)
            } else {
                examples <- character(0)
            }
            examples <- examples[!is.na(examples) & nzchar(examples)]
            if (length(examples) > 0) {
                drop.examples[i] <- paste(head(examples, 3), collapse = " | ")
            }
        }
    }

    out <- data.frame(time = df$t,
                      d = df$d,
                      senspar = df$sens.par,
                      n.drop.sets = n.drop.sets,
                      drop.examples = drop.examples)
    rownames(out) <- NULL
    out
}

.npsa_rv_summary_table <- function(res.RV) {
    if (is.null(res.RV) || is.null(res.RV$res.table)) return(NULL)
    df <- as.data.frame(res.RV$res.table)
    if (nrow(df) == 0) return(NULL)

    out <- data.frame(time = df$t0,
                      theta = df$theta,
                      rho = if ("rho" %in% names(df)) df$rho else NA,
                      RV = df$RV,
                      MIRV = df$MIRV,
                      sp.RV = .npsa_sp_from_rv(df$RV),
                      sp.MIRV = .npsa_sp_from_rv(df$MIRV),
                      lower.b = if ("lower.b" %in% names(df)) df$lower.b else NA)
    rownames(out) <- NULL
    out
}

.npsa_urv_summary_table <- function(res.RV, fit.times = NULL) {
    if (is.null(res.RV) || is.null(res.RV$unif.RV)) return(NULL)

    if (!is.null(res.RV$uniform.window)) {
        t.lower <- res.RV$uniform.window[1]
        t.upper <- res.RV$uniform.window[2]
    } else if (!is.null(res.RV$unif.idx) && !is.null(fit.times)) {
        t.lower <- min(fit.times[res.RV$unif.idx], na.rm = TRUE)
        t.upper <- max(fit.times[res.RV$unif.idx], na.rm = TRUE)
    } else {
        t.lower <- NA
        t.upper <- NA
    }

    out <- data.frame(t.lower = t.lower,
                      t.upper = t.upper,
                      URV = res.RV$unif.RV,
                      sp.URV = .npsa_sp_from_rv(res.RV$unif.RV),
                      window.source = if (is.null(res.RV$uniform.window.source)) NA else res.RV$uniform.window.source)
    rownames(out) <- NULL
    out
}

.npsa_time_summary_table <- function(object) {
    time.info <- object$time.info
    if (is.null(time.info)) return(NULL)

    show.times <- function(x) {
        if (is.null(x)) return(NA_character_)
        if (length(x) == 0) return("")
        x <- signif(x, 4)
        if (length(x) <= 6) {
            paste(x, collapse = ", ")
        } else {
            paste0(paste(x[1:6], collapse = ", "), ", ...")
        }
    }

    out <- data.frame(
        item = c("fit.times", "report.times", "rv.times", "uniform.window"),
        n = c(length(time.info$fit.times),
              length(time.info$report.times),
              length(time.info$rv.times),
              length(time.info$uniform.window)),
        values = c(show.times(time.info$fit.times),
                   show.times(time.info$report.times),
                   show.times(time.info$rv.times),
                   show.times(time.info$uniform.window))
    )
    rownames(out) <- NULL
    out
}

.npsa_surv_summary_tables <- function(object) {
    report.times <- object$time.info$report.times
    out <- list(time.settings = .npsa_time_summary_table(object))

    if (!is.null(object$bounds.df) && !is.null(object$bounds.df$bounds.df)) {
        bounds.raw <- object$bounds.df$bounds.df
        out$surv.diff <- .npsa_bounds_summary_table(
            bounds.raw[bounds.raw$d == 0, , drop = FALSE],
            report.times = report.times,
            effect.name = "surv.diff",
            include.d = FALSE
        )
        out$bounds <- .npsa_bounds_summary_table(
            bounds.raw[bounds.raw$d != 0, , drop = FALSE],
            report.times = report.times,
            effect.name = "surv.diff",
            include.d = TRUE
        )
    }

    out$senspar <- .npsa_senspar_summary_table(object$senspar.df, report.times)
    out$rv <- .npsa_rv_summary_table(object$res.RV)
    out$urv <- .npsa_urv_summary_table(object$res.RV, object$result$fit.times)

    if (!is.null(object$bounds.df) && !is.null(object$bounds.df$bounds.df.rmst)) {
        out$rmst <- .npsa_bounds_summary_table(
            object$bounds.df$bounds.df.rmst,
            report.times = NULL,
            effect.name = "rmst.diff",
            include.d = TRUE
        )
    }

    out
}

.npsa_print_summary_table <- function(title, tbl, digits = 3,
                                      max.rows = 12, object.name = NULL) {
    cat("\n", title, ":\n", sep = "")
    if (is.null(tbl) || nrow(tbl) == 0) {
        cat("Not available.\n")
        return(invisible(NULL))
    }

    tbl.print <- tbl
    n.more <- 0
    if (is.finite(max.rows) && nrow(tbl.print) > max.rows) {
        n.more <- nrow(tbl.print) - max.rows
        tbl.print <- tbl.print[seq_len(max.rows), , drop = FALSE]
    }

    num.cols <- sapply(tbl.print, is.numeric)
    tbl.print[, num.cols] <- lapply(tbl.print[, num.cols, drop = FALSE], function(x) round(x, digits))
    print(tbl.print, row.names = FALSE)

    if (n.more > 0) {
        if (is.null(object.name)) {
            cat("... ", n.more, " more rows not printed.\n", sep = "")
        } else {
            cat("... ", n.more, " more rows not printed. See ", object.name, ".\n", sep = "")
        }
    }
    invisible(tbl)
}

#' Summarize NPSA Survival Sensitivity Results
#'
#' @param object An object returned by \code{\link{npsa_surv}()}.
#' @param type Which summary to print. Options are \code{"overview"},
#'   \code{"bounds"}, \code{"senspar"}, \code{"rv"}, \code{"rmst"}, and
#'   \code{"all"}.
#' @param digits Number of digits for printing.
#' @param ... Additional arguments.
#'
#' @return Invisibly returns \code{object}. Clean user-facing tables are also
#'   stored in \code{object$summary.tables}.
#'
#' @export
#' @method summary npsa_surv
summary.npsa_surv <- function(object,
                              type = c("overview", "bounds", "senspar", "rv", "rmst", "all"),
                              digits = 3, ...) {
    type <- match.arg(type)

    tables <- object$summary.tables
    if (is.null(tables)) tables <- .npsa_surv_summary_tables(object)

    cat("NPSA Survival Sensitivity Report\n")
    cat("--------------------------------\n")

    max.rows <- if (type == "all") Inf else 12

    if (type %in% c("overview", "all")) {
        .npsa_print_summary_table("Time settings", tables$time.settings, digits, max.rows,
                                  "object$summary.tables$time.settings")
        .npsa_print_summary_table("No-unobserved-confounding survival difference",
                                  tables$surv.diff, digits, max.rows,
                                  "object$summary.tables$surv.diff")

        if (!is.null(object$senspar.df$meta)) {
            cat("\nSensitivity parameter simulation:\n")
            cat("rep:", object$senspar.df$meta$rep, "\n")
            cat("seed:", object$senspar.df$meta$seed, "\n")
            if (!is.null(tables$senspar)) {
                cat("d:", paste(unique(tables$senspar$d), collapse = ", "), "\n")
            }
        }
        .npsa_print_summary_table("Sensitivity parameter summary", tables$senspar,
                                  digits, max.rows, "object$summary.tables$senspar")
        .npsa_print_summary_table("Sensitivity bounds", tables$bounds,
                                  digits, max.rows, "object$summary.tables$bounds")
        .npsa_print_summary_table("RV and MIRV", tables$rv,
                                  digits, max.rows, "object$summary.tables$rv")
        .npsa_print_summary_table("Uniform RV", tables$urv,
                                  digits, max.rows, "object$summary.tables$urv")
        if (!is.null(tables$rmst)) {
            .npsa_print_summary_table("RMST difference", tables$rmst,
                                      digits, max.rows, "object$summary.tables$rmst")
        }
    } else if (type == "bounds") {
        .npsa_print_summary_table("No-unobserved-confounding survival difference",
                                  tables$surv.diff, digits, max.rows,
                                  "object$summary.tables$surv.diff")
        .npsa_print_summary_table("Sensitivity bounds", tables$bounds,
                                  digits, max.rows, "object$summary.tables$bounds")
    } else if (type == "senspar") {
        .npsa_print_summary_table("Sensitivity parameter summary", tables$senspar,
                                  digits, max.rows, "object$summary.tables$senspar")
    } else if (type == "rv") {
        .npsa_print_summary_table("RV and MIRV", tables$rv,
                                  digits, max.rows, "object$summary.tables$rv")
        .npsa_print_summary_table("Uniform RV", tables$urv,
                                  digits, max.rows, "object$summary.tables$urv")
    } else if (type == "rmst") {
        .npsa_print_summary_table("RMST difference", tables$rmst,
                                  digits, max.rows, "object$summary.tables$rmst")
    }

    invisible(object)
}



#' Summarize No-Unobserved-Confounding Survival Results
#'
#' @param object An object returned by \code{\link{np_surv}()}.
#' @param digits Number of digits for printing.
#' @param ... Additional arguments.
#'
#' @return Invisibly returns \code{object}.
#' If \code{np_surv(rmst = TRUE)} was used, the summary also prints the RMST
#' difference table.
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

    if (!is.null(object$rmst.summary)) {
        cat("\nRMST difference summary:\n")
        tbl <- object$rmst.summary
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
#'   \code{"surv.ratio"}, \code{"risk.ratio"}, and \code{"nnt"}.
#' @param uniform Logical; whether to add uniform confidence bands.
#' @param transform Logical; for treatment-specific survival curves, whether to
#'   use logit-scale pointwise confidence intervals and logit/equi-precision
#'   uniform bands. Contrast plots use their CFsurvival-style inference scale.
#' @param ... Additional arguments.
#'
#' @return A \code{ggplot} object.
#'
#' @export
#' @method plot npSurv
plot.npSurv <- function(x, type = c("surv", "surv.diff", "surv.ratio", "risk.ratio", "nnt"),
                        uniform = TRUE, transform = TRUE, ...) {
    type <- match.arg(type)

    if (type == "surv") {
        return(.plot.np_surv_curves(x$surv.df, uniform = uniform, transform = transform))
    }
    return(.plot.np_contrast(x, type = type, uniform = uniform))
}

.plot.np_surv_curves <- function(df, uniform = TRUE, transform = TRUE) {
    if (transform) {
        ptwise.lower <- "ptwise.logit.lower"
        ptwise.upper <- "ptwise.logit.upper"
        unif.lower <- "unif.logit.lower"
        unif.upper <- "unif.logit.upper"
    } else {
        ptwise.lower <- "ptwise.lower"
        ptwise.upper <- "ptwise.upper"
        unif.lower <- "unif.ew.lower"
        unif.upper <- "unif.ew.upper"
    }

    needed <- c(ptwise.lower, ptwise.upper)
    if (!all(needed %in% names(df))) {
        stop("Requested pointwise confidence interval columns were not found.")
    }

    df$ptwise.plot.lower <- df[[ptwise.lower]]
    df$ptwise.plot.upper <- df[[ptwise.upper]]

    has.uniform <- isTRUE(uniform) && all(c(unif.lower, unif.upper) %in% names(df))
    if (has.uniform) {
        df$uniform.plot.lower <- df[[unif.lower]]
        df$uniform.plot.upper <- df[[unif.upper]]
    }

    p <- ggplot(df, aes(x = time, y = surv, color = as.factor(trt), group = trt)) +
        geom_step() +
        geom_step(aes(y = ptwise.plot.lower), linetype = "dashed", na.rm = TRUE) +
        geom_step(aes(y = ptwise.plot.upper), linetype = "dashed", na.rm = TRUE) +
        scale_color_manual(values = c("0" = "#0072B2", "1" = "#D55E00"),
                           labels = c("0" = "Control", "1" = "Treatment")) +
        labs(color = "Treatment") +
        xlab("Time") +
        ylab("Treatment-specific survival") +
        coord_cartesian(ylim = c(0, 1)) +
        theme_bw() +
        theme(legend.position = "bottom",
              legend.title = element_blank(),
              panel.grid.minor = element_blank())

    if (has.uniform) {
        p <- p +
            geom_step(aes(y = uniform.plot.lower), linetype = "longdash", na.rm = TRUE) +
            geom_step(aes(y = uniform.plot.upper), linetype = "longdash", na.rm = TRUE)
    }
    return(p)
}

.plot.np_contrast <- function(x, type, uniform = TRUE) {
    if (type == "surv.diff" && "times" %in% names(x$surv.diff.df)) {
        return(.plot.np_surv_diff_bounds(x$surv.diff.df, uniform = uniform))
    }

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
    df$estimate <- df[[info$est]]

    p <- ggplot(df, aes(x = time)) +
        geom_line(aes(y = estimate, color = "Estimate"), na.rm = TRUE) +
        geom_line(aes(y = ptwise.lower, color = "Pointwise CI"), linetype = "dashed", na.rm = TRUE) +
        geom_line(aes(y = ptwise.upper, color = "Pointwise CI"), linetype = "dashed", na.rm = TRUE) +
        scale_color_manual(values = c("Estimate" = "black",
                                      "Pointwise CI" = "#0072B2",
                                      "Uniform Band" = "#009E73")) +
        xlab("Time") +
        ylab(info$ylab) +
        theme_bw() +
        theme(legend.position = "bottom",
              legend.title = element_blank(),
              panel.grid.minor = element_blank())

    if (isTRUE(uniform)) {
        p <- p +
            geom_line(aes(y = unif.lower, color = "Uniform Band"), linetype = "longdash", na.rm = TRUE) +
            geom_line(aes(y = unif.upper, color = "Uniform Band"), linetype = "longdash", na.rm = TRUE)
    }

    if (is.finite(info$ref)) {
        p <- p + geom_hline(yintercept = info$ref, color = "grey45", linetype = "dotted")
    }
    return(p)
}

.plot.np_surv_diff_bounds <- function(df, uniform = TRUE) {
    if ("ptwise.trans.lower" %in% names(df)) {
        ptwise.lower <- "ptwise.trans.lower"
        ptwise.upper <- "ptwise.trans.upper"
        unif.lower <- "uniform.trans.lower"
        unif.upper <- "uniform.trans.upper"
    } else {
        ptwise.lower <- "ptwise.bounds.lower"
        ptwise.upper <- "ptwise.bounds.upper"
        unif.lower <- "uniform.bounds.lower"
        unif.upper <- "uniform.bounds.upper"
    }

    df$ptwise.plot.lower <- df[[ptwise.lower]]
    df$ptwise.plot.upper <- df[[ptwise.upper]]
    df$uniform.plot.lower <- df[[unif.lower]]
    df$uniform.plot.upper <- df[[unif.upper]]

    p <- ggplot(df, aes(x = times)) +
        geom_line(aes(y = theta.obs, color = "Estimate"), na.rm = TRUE) +
        geom_line(aes(y = ptwise.plot.lower, color = "Pointwise CI"),
                  linetype = "dashed", na.rm = TRUE) +
        geom_line(aes(y = ptwise.plot.upper, color = "Pointwise CI"),
                  linetype = "dashed", na.rm = TRUE) +
        geom_hline(yintercept = 0, color = "grey45", linetype = "dotted") +
        scale_color_manual(values = c("Estimate" = "black",
                                      "Pointwise CI" = "#0072B2",
                                      "Uniform Band" = "#009E73")) +
        xlab("Time") +
        ylab("Survival difference (treatment - control)") +
        theme_bw() +
        theme(legend.position = "bottom",
              legend.title = element_blank(),
              panel.grid.minor = element_blank())

    if (isTRUE(uniform)) {
        p <- p +
            geom_line(aes(y = uniform.plot.lower, color = "Uniform Band"),
                      linetype = "longdash", na.rm = TRUE) +
            geom_line(aes(y = uniform.plot.upper, color = "Uniform Band"),
                      linetype = "longdash", na.rm = TRUE)
    }
    return(p)
}

#' Estimate observed components
#'
#' @keywords internal
.report.bounds <- function(bound.times, result, rho=1, band.end.pts = c(0,Inf), conf.level=.95, boot=10000,
                           sens.df.mean = NULL, num_drop = NULL, pct_drop = NULL, n_var = NULL,
                           rmst = TRUE, sens.rmst.df.mean = NULL, transform = TRUE, scale = TRUE) {

    if (is.null(bound.times)) bound.times <- result$fit.times
    if (any(bound.times > max(result$fit.times))) {
        message("Some bound.times > maximum observed event time - removed for bounds.")
        bound.times <- bound.times[bound.times <= max(result$fit.times)]
    }

    if (is.null(sens.df.mean)) {

        obs.est.idx <- sapply(bound.times, function(x) {
            which(near(x, result$fit.times))
        })

        # Estimate lower and upper bounds without sensitivity (assume zero sensitivity)
        effect.bounds <- .get.effect.bounds(
            fit.times = result$fit.times[obs.est.idx],
            theta.obs = result$obs.comps.df$theta.obs[obs.est.idx],
            psi = result$obs.comps.df$psi[obs.est.idx],
            tau = result$tau,
            sens.out = rep(0, length(bound.times)),
            sens.trt = 0,
            rho = rho
        )

        bounds.conf.int <- .bounds.confints(
            effect.bounds,
            psi = result$obs.comps.df$psi[obs.est.idx],
            tau = result$tau,
            IF.vals.theta.obs = result$IF.vals.theta.obs[, obs.est.idx, drop = FALSE],
            IF.vals.psi = result$IF.vals.psi[, obs.est.idx, drop = FALSE],
            IF.vals.tau = result$IF.vals.tau,
            rho = rho,
            band.end.pts = band.end.pts,
            conf.level = conf.level,
            scale = scale,
            boot = boot
        )

        bounds.df <- bounds2df(bounds.conf.int, theta.obs = result$obs.comps.df$theta.obs[obs.est.idx],
                               d = NULL, transform = transform, time.zero = TRUE,
                               effect.range = c(-1, 1))

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

            senspar.idx <- sapply(bound.times, function(x) {
                which(near(x, sens.df.mean$t[sens.df.mean$d == d]))
            })

            obs.est.idx <- sapply(bound.times, function(x) {
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
                IF.vals.theta.obs = result$IF.vals.theta.obs[, obs.est.idx, drop = FALSE],
                IF.vals.psi = result$IF.vals.psi[, obs.est.idx, drop = FALSE],
                IF.vals.tau = result$IF.vals.tau,
                conf.level = conf.level,
                scale = scale
            )

            df <- bounds2df(bounds.conf.int, theta.obs = result$obs.comps.df$theta.obs[obs.est.idx],
                            d = d, transform = transform, time.zero = TRUE,
                            effect.range = c(-1, 1))

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

#' Plot NPSA Survival Sensitivity Bounds
#'
#' Plot method for objects returned by \code{\link{npsa_surv}()}.
#'
#' @param x An object returned by \code{\link{npsa_surv}()}.
#' @param ... Additional arguments passed to \code{\link{plot.boundsdf}()}.
#'
#' @return A \code{ggplot} object showing sensitivity bounds.
#'
#' @export
#' @method plot npsa_surv
plot.npsa_surv <- function(x, ...) {
    if (is.null(x$bounds.df) || is.null(x$bounds.df$bounds.df)) {
        stop("`x` does not contain bounds results. If `senspar.only = TRUE`, run `npsa_surv()` again without `senspar.only` before plotting.")
    }
    plot(x$bounds.df$bounds.df, ...)
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
