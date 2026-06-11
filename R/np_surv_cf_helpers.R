# CFsurvival helper logic adapted from tedwestling/CFsurvival
# R/fit_survival.R and R/get_comparisons.R.
# See inst/CFSurvival_vendor_notes.txt for provenance and local changes.

.np_estimate.uniform.quantile <- function(IF.vals, conf.level=.95, scale=TRUE) {
    IF.vals <- cbind(IF.vals)
    n <- nrow(IF.vals)
    if(scale) IF.vals <- scale(IF.vals)
    maxes <- replicate(1e4, max(abs(rbind(stats::rt(n, df = n - 1)/sqrt(n)) %*% IF.vals)))
    return(list(quantile=stats::quantile(maxes, conf.level), maxes=maxes))
}

.np_surv.confints <- function(times, est, IF.vals, isotonize=TRUE, conf.band=TRUE,
                              band.end.pts=c(0,Inf), conf.level=.95,
                              band.ew=FALSE) {
    logit <- function(x) log(x / (1-x))
    logit.prime <- function(x) 1/(x * (1-x))
    expit <- function(x) 1/(1 + exp(-x))

    n <- nrow(IF.vals)
    IF.vals.logit <- IF.vals
    for(j in 1:length(est)) IF.vals.logit[,j] <- IF.vals[,j] * logit.prime(est[j])

    res <- NULL
    res$se <- sqrt(colMeans(IF.vals^2)) / sqrt(n)
    res$se.logit <- sqrt(colMeans(IF.vals.logit^2)) / sqrt(n)

    res$se[res$se == 0] <- NA
    res$se.logit[is.infinite(res$se.logit)] <- NA
    quant <- stats::qt(1-(1-conf.level)/2, n - 1)

    res$ptwise.lower <- pmax(est - quant * res$se, 0)
    res$ptwise.upper <- pmin(est + quant * res$se, 1)
    res$ptwise.logit.lower <- expit(logit(est) - quant * res$se.logit)
    res$ptwise.logit.upper <- expit(logit(est) + quant * res$se.logit)

    out <- NULL
    if(conf.band) {
        ew.idx <- !is.na(res$se)
        if (band.ew) ew.idx <- ew.idx & times >= band.end.pts[1] & times <= band.end.pts[2]
        if(any(ew.idx)) {
            unif.vals <- .np_estimate.uniform.quantile(IF.vals[,ew.idx, drop=FALSE],
                                                       conf.level, scale = FALSE)
            unif.quant <- unif.vals$quantile
            out$ew.sim.maxes <- unif.vals$maxes
            out$unif.ew.quant <- unif.quant
            res$unif.ew.lower <- pmax(est - unif.quant / sqrt(n), 0)
            res$unif.ew.upper <- pmin(est + unif.quant / sqrt(n), 1)
            if (band.ew) {
                res$unif.ew.lower[times < band.end.pts[1] | times > band.end.pts[2]] <- NA
                res$unif.ew.upper[times < band.end.pts[1] | times > band.end.pts[2]] <- NA
            }

            if(isotonize) {
                res$unif.ew.lower[!is.na(res$unif.ew.lower)] <-
                    1 - stats::isoreg(times[!is.na(res$unif.ew.lower)],
                                      1-res$unif.ew.lower[!is.na(res$unif.ew.lower)])$yf
                res$unif.ew.upper[!is.na(res$unif.ew.upper)] <-
                    1 - stats::isoreg(times[!is.na(res$unif.ew.upper)],
                                      1-res$unif.ew.upper[!is.na(res$unif.ew.upper)])$yf
            }
        } else {
            res$unif.ew.lower <- rep(NA, length(res$se))
            res$unif.ew.upper <- rep(NA, length(res$se))
        }

        if(any(!is.na(res$se.logit) & times >= band.end.pts[1] & times <= band.end.pts[2]))  {
            unif.logit.vals <-
                .np_estimate.uniform.quantile(IF.vals.logit[,!is.na(res$se.logit) &
                                                                times >= band.end.pts[1] &
                                                                times <= band.end.pts[2],
                                                            drop=FALSE],
                                             conf.level)
            unif.logit.quant <- unif.logit.vals$quantile
            out$logit.sim.maxes <- unif.logit.vals$maxes
            out$unif.logit.quant <- unif.logit.quant
            res$unif.logit.lower <- expit(logit(est) - unif.logit.quant * res$se.logit)
            res$unif.logit.upper <- expit(logit(est) + unif.logit.quant * res$se.logit)
            res$unif.logit.lower[times < band.end.pts[1] | times > band.end.pts[2]] <- NA
            res$unif.logit.upper[times < band.end.pts[1] | times > band.end.pts[2]] <- NA

            if(isotonize) {
                res$unif.logit.lower[!is.na(res$unif.logit.lower)] <-
                    1 - stats::isoreg(times[!is.na(res$unif.logit.lower)],
                                      1-res$unif.logit.lower[!is.na(res$unif.logit.lower)])$yf
                res$unif.logit.upper[!is.na(res$unif.logit.upper)] <-
                    1 - stats::isoreg(times[!is.na(res$unif.logit.upper)],
                                      1-res$unif.logit.upper[!is.na(res$unif.logit.upper)])$yf
            }
        } else {
            res$unif.logit.lower <-  rep(NA, length(res$se))
            res$unif.logit.upper <-  rep(NA, length(res$se))
        }
    } else {
        res$unif.ew.lower <- res$unif.ew.upper <- rep(NA, length(res$se))
        res$unif.logit.lower <- res$unif.logit.upper <- rep(NA, length(res$se))
    }
    out$res <- res
    return(out)
}

.np_surv.difference <- function(fit.times, surv.0, surv.1, IF.vals.0, IF.vals.1,
                                conf.band=TRUE, band.end.pts=c(0,Inf), conf.level=.95) {
    logit <- function(x) log(x / (1-x))
    logit.prime <- function(x) 1/(x * (1-x))
    expit <- function(x) 1/(1 + exp(-x))

    n <- nrow(IF.vals.0)
    quant <- stats::qt(1-(1-conf.level)/2, n - 1)
    df <- data.frame(time=c(0,fit.times), surv.diff=c(0,surv.1-surv.0))

    diff <- surv.1 - surv.0
    IF.diff <- IF.vals.1 - IF.vals.0
    se.diff <- sqrt(colMeans(IF.diff^2))
    se.diff[se.diff == 0] <- NA

    logit.diff <- logit((diff + 1)/2)
    IF.diff.logit <- IF.diff
    for(j in 1:length(diff)) IF.diff.logit[,j] <- IF.diff[,j] * logit.prime((diff[j] + 1) / 2) / 2

    se.diff.logit <- sqrt(colMeans(IF.diff.logit^2))
    se.diff.logit[se.diff.logit == 0] <- NA

    df$se <- c(0,se.diff)/sqrt(n)
    df$se.logit <- c(0, se.diff.logit / sqrt(n))
    ll <- 2 * expit(logit.diff - quant * se.diff.logit / sqrt(n)) - 1
    ul <- 2 * expit(logit.diff + quant * se.diff.logit / sqrt(n)) - 1
    df$ptwise.lower <- c(0, ll)
    df$ptwise.upper <- c(0, ul)
    df$ptwise.pval <- c(1, stats::pchisq((logit.diff / (se.diff.logit / sqrt(n)))^2,
                                         df=1, lower.tail = FALSE))

    has.band <- conf.band & any(!is.na(se.diff.logit) &
                                fit.times >= band.end.pts[1] & fit.times <= band.end.pts[2])
    if(has.band) {
        unif.info <- .np_estimate.uniform.quantile(IF.diff.logit[,!is.na(se.diff.logit) &
                                                                    fit.times >= band.end.pts[1] &
                                                                    fit.times <= band.end.pts[2],
                                                                drop=FALSE],
                                                  conf.level)
        unif.quant <- unif.info$quantile

        ll <- 2 * expit(logit.diff - unif.quant * se.diff.logit / sqrt(n)) - 1
        ul <- 2 * expit(logit.diff + unif.quant * se.diff.logit / sqrt(n)) - 1
        df$unif.lower <- c(0, ll)
        df$unif.upper <- c(0, ul)

        df$unif.lower[df$time < band.end.pts[1] | df$time > band.end.pts[2]] <- NA
        df$unif.upper[df$time < band.end.pts[1] | df$time > band.end.pts[2]] <- NA
    } else {
        df$unif.lower <- df$unif.upper <- NA
    }

    res <- list(surv.diff.df=df)
    if(has.band) {
        res$surv.diff.unif.quant <- unif.quant
        res$surv.diff.sim.maxes <- unif.info$maxes
    }
    return(res)
}

.np_surv.ratio <- function(fit.times, surv.0, surv.1, IF.vals.0, IF.vals.1,
                           conf.band=TRUE, band.end.pts=c(0, Inf),conf.level=.95) {
    n <- nrow(IF.vals.0)
    quant <- stats::qt(1-(1-conf.level)/2, df = n - 1)
    df <- data.frame(time=c(0,fit.times), log.surv.ratio=c(0,log(surv.1) - log(surv.0)),
                     surv.ratio=c(1,surv.1 / surv.0))

    IF.log.ratio <- IF.vals.1 / matrix(surv.1, nrow=nrow(IF.vals.1), ncol=ncol(IF.vals.1), byrow=TRUE) -
        IF.vals.0 / matrix(surv.0, nrow=nrow(IF.vals.0), ncol=ncol(IF.vals.0), byrow=TRUE)
    se.log.ratio <- sqrt(colMeans(IF.log.ratio^2))
    se.log.ratio[se.log.ratio == 0] <- NA
    se.log.ratio[is.infinite(se.log.ratio)] <- NA
    df$se.log.ratio <- c(0,se.log.ratio)/sqrt(n)
    df$ptwise.lower.log <- df$log.surv.ratio - quant * df$se.log.ratio
    df$ptwise.upper.log <- df$log.surv.ratio + quant * df$se.log.ratio
    df$ptwise.lower <- exp(df$ptwise.lower.log)
    df$ptwise.upper <- exp(df$ptwise.upper.log)
    df$ptwise.pval <- c(1,stats::pchisq((df$log.surv.ratio[-1]/(df$se.log.ratio[-1]))^2,
                                        df=1, lower.tail = FALSE))

    has.band <- conf.band & any(!is.na(se.log.ratio) &
                                fit.times >= band.end.pts[1] & fit.times <= band.end.pts[2])
    if(has.band) {
        log.unif.info <- .np_estimate.uniform.quantile(IF.log.ratio[,!is.na(se.log.ratio) &
                                                                       fit.times >= band.end.pts[1] &
                                                                       fit.times <= band.end.pts[2],
                                                                   drop=FALSE],
                                                      conf.level)
        log.unif.quant <- log.unif.info$quantile
        df$unif.lower.log <- df$log.surv.ratio - log.unif.quant * df$se.log.ratio
        df$unif.upper.log <- df$log.surv.ratio + log.unif.quant * df$se.log.ratio
        df$unif.lower <- exp(df$unif.lower.log)
        df$unif.upper <- exp(df$unif.upper.log)

        df$unif.lower[df$time < band.end.pts[1] | df$time > band.end.pts[2]] <- NA
        df$unif.upper[df$time < band.end.pts[1] | df$time > band.end.pts[2]] <- NA
    } else {
        df$unif.lower.log <- df$unif.upper.log <- df$unif.lower <- df$unif.upper <- NA
    }

    res <- list(surv.ratio.df=df)
    if(has.band)  {
        res$log.surv.ratio.unif.quant <- log.unif.quant
        res$log.surv.ratio.sim.maxes <- log.unif.info$maxes
    }

    return(res)
}

.np_risk.ratio <- function(fit.times, surv.0, surv.1, IF.vals.0, IF.vals.1,
                           conf.band=TRUE, band.end.pts=c(0,Inf), conf.level=.95) {
    n <- nrow(IF.vals.0)
    quant <- stats::qt(1-(1-conf.level)/2, df = n - 1)
    risk.0 <- 1-surv.0
    risk.1 <- 1-surv.1
    IF.vals.0 <- -IF.vals.0
    IF.vals.1 <- -IF.vals.1
    df <- data.frame(time=c(0,fit.times), log.risk.ratio=c(NA,log(risk.1) - log(risk.0)),
                     risk.ratio=c(NA,risk.1 / risk.0))

    IF.log.ratio <- IF.vals.1 / matrix(risk.1, nrow=nrow(IF.vals.1), ncol=ncol(IF.vals.1), byrow=TRUE) -
        IF.vals.0 / matrix(risk.0, nrow=nrow(IF.vals.0), ncol=ncol(IF.vals.0), byrow=TRUE)
    se.log.ratio <- sqrt(colMeans(IF.log.ratio^2))
    se.log.ratio[se.log.ratio == 0] <- NA
    se.log.ratio[is.infinite(se.log.ratio)] <- NA
    df$se.log.ratio <- c(NA,se.log.ratio)/sqrt(n)
    df$ptwise.lower.log <- df$log.risk.ratio - quant * df$se.log.ratio
    df$ptwise.upper.log <- df$log.risk.ratio + quant * df$se.log.ratio
    df$ptwise.lower <- exp(df$ptwise.lower.log)
    df$ptwise.upper <- exp(df$ptwise.upper.log)
    df$ptwise.pval <- c(1,stats::pchisq((df$log.risk.ratio[-1]/(df$se.log.ratio[-1]))^2,
                                        df=1, lower.tail = FALSE))

    has.band <- conf.band & any(!is.na(se.log.ratio) &
                                fit.times >= band.end.pts[1] & fit.times <= band.end.pts[2])
    if(has.band) {
        log.unif.info <- .np_estimate.uniform.quantile(IF.log.ratio[,!is.na(se.log.ratio) &
                                                                       fit.times >= band.end.pts[1] &
                                                                       fit.times <= band.end.pts[2],
                                                                   drop=FALSE],
                                                      conf.level)
        log.unif.quant <- log.unif.info$quantile
        df$unif.lower.log <- df$log.risk.ratio - log.unif.quant * df$se.log.ratio
        df$unif.upper.log <- df$log.risk.ratio + log.unif.quant * df$se.log.ratio
        df$unif.lower <- exp(df$unif.lower.log)
        df$unif.upper <- exp(df$unif.upper.log)

        df$unif.lower[df$time < band.end.pts[1] | df$time > band.end.pts[2]] <- NA
        df$unif.upper[df$time < band.end.pts[1] | df$time > band.end.pts[2]] <- NA
    } else {
        df$unif.lower.log <- df$unif.upper.log <- df$unif.lower <- df$unif.upper <- NA
    }

    res <- list(risk.ratio.df=df)
    if(has.band) {
        res$log.risk.ratio.unif.quant <- log.unif.quant
        res$log.risk.ratio.sim.maxes <- log.unif.info$maxes
    }

    return(res)
}

.np_nnt <- function(fit.times, surv.0, surv.1, IF.vals.0, IF.vals.1,
                    conf.band=TRUE, band.end.pts=c(0, Inf), conf.level=.95) {
    n <- nrow(IF.vals.0)
    quant <- stats::qt(1-(1-conf.level)/2, df = n - 1)
    df <- data.frame(time=c(0,fit.times), nnt=c(NA,1/(surv.1-surv.0)))

    IF.nnt <- (IF.vals.0 - IF.vals.1) /
        matrix((surv.1 - surv.0)^2, nrow=nrow(IF.vals.0), ncol=ncol(IF.vals.0), byrow=TRUE)
    se.nnt <- sqrt(colMeans(IF.nnt^2))
    se.nnt[se.nnt == 0] <- NA

    df$se <- c(0,se.nnt)/sqrt(n)
    df$ptwise.lower <- df$nnt - quant * df$se
    df$ptwise.upper <- df$nnt + quant * df$se

    has.band <- conf.band & any(!is.na(se.nnt) & !is.nan(se.nnt) & !is.infinite(se.nnt) &
                                fit.times >= band.end.pts[1] & fit.times <= band.end.pts[2])
    if(has.band) {
        unif.info <- .np_estimate.uniform.quantile(IF.nnt[,!is.na(se.nnt) &
                                                             !is.nan(se.nnt) &
                                                             !is.infinite(se.nnt) &
                                                             fit.times >= band.end.pts[1] &
                                                             fit.times <= band.end.pts[2],
                                                         drop=FALSE],
                                                  conf.level)
        unif.quant <- unif.info$quantile
        df$unif.lower <- df$nnt - unif.quant * df$se
        df$unif.upper <- df$nnt + unif.quant * df$se
        df$unif.lower[df$time < band.end.pts[1] | df$time > band.end.pts[2]] <- NA
        df$unif.upper[df$time < band.end.pts[1] | df$time > band.end.pts[2]] <- NA
    } else {
        df$unif.lower <- df$unif.upper <- NA
    }

    res <- list(nnt.df=df)
    if(has.band) {
        res$nnt.unif.quant <- unif.quant
        res$nnt.sim.maxes <- unif.info$maxes
    }
    return(res)
}

.np_report_cf_surv <- function(time, event, treat, result, conf.band=TRUE,
                               conf.level=.95, contrasts=c("surv.diff", "surv.ratio"),
                               uniform.cutpoint=c(0.01, 0.99), uniform.window=NULL,
                               isotonize=TRUE) {
    fit.times <- result$fit.times
    surv.0 <- .np_get_surv_object(result, trt = 0, isotonize = isotonize)
    surv.1 <- .np_get_surv_object(result, trt = 1, isotonize = isotonize)
    exact.window <- !is.null(uniform.window)

    if (is.null(uniform.window)) {
        band.end.pts.0 <- .np_surv_band_endpts(time[event == 1 & treat == 0],
                                               surv.0$surv.iso, fit.times,
                                               uniform.cutpoint)
        band.end.pts.1 <- .np_surv_band_endpts(time[event == 1 & treat == 1],
                                               surv.1$surv.iso, fit.times,
                                               uniform.cutpoint)
        band.end.pts <- .np_contrast_band_endpts(time[event == 1],
                                                 surv.0$surv.iso, surv.1$surv.iso,
                                                 fit.times, uniform.cutpoint)
    } else {
        band.end.pts <- .np_trim_uniform_window(uniform.window, fit.times)
        band.end.pts.0 <- band.end.pts
        band.end.pts.1 <- band.end.pts
    }

    surv.df.0 <- .np_surv_df_one(fit.times, surv.0, trt = 0, conf.band = conf.band,
                                 band.end.pts = band.end.pts.0,
                                 conf.level = conf.level,
                                 isotonize = isotonize,
                                 band.ew = exact.window)
    surv.df.1 <- .np_surv_df_one(fit.times, surv.1, trt = 1, conf.band = conf.band,
                                 band.end.pts = band.end.pts.1,
                                 conf.level = conf.level,
                                 isotonize = isotonize,
                                 band.ew = exact.window)

    out <- list(surv.df = rbind(surv.df.0$surv.df, surv.df.1$surv.df),
                surv.0.unif.ew.quant = surv.df.0$unif.ew.quant,
                surv.0.unif.logit.quant = surv.df.0$unif.logit.quant,
                surv.0.ew.sim.maxes = surv.df.0$ew.sim.maxes,
                surv.0.logit.sim.maxes = surv.df.0$logit.sim.maxes,
                surv.1.unif.ew.quant = surv.df.1$unif.ew.quant,
                surv.1.unif.logit.quant = surv.df.1$unif.logit.quant,
                surv.1.ew.sim.maxes = surv.df.1$ew.sim.maxes,
                surv.1.logit.sim.maxes = surv.df.1$logit.sim.maxes)

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
                            band.end.pts=c(0, Inf), conf.level=.95,
                            isotonize=TRUE, band.ew=FALSE) {
    c.int <- .np_surv.confints(fit.times, surv$surv, surv$IF.vals,
                               conf.band = conf.band,
                               band.end.pts = band.end.pts,
                               conf.level = conf.level,
                               isotonize = isotonize,
                               band.ew = band.ew)
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
                unif.logit.quant = c.int$unif.logit.quant,
                ew.sim.maxes = c.int$ew.sim.maxes,
                logit.sim.maxes = c.int$logit.sim.maxes))
}

.np_surv_band_endpts <- function(event.times, surv.iso, fit.times, uniform.cutpoint) {
    event.times <- event.times[is.finite(event.times)]
    if (length(event.times) == 0) return(c(min(fit.times), max(fit.times)))

    lower <- as.numeric(stats::quantile(event.times, uniform.cutpoint[1], na.rm = TRUE))
    surv.cut <- 1 - uniform.cutpoint[2]
    upper.idx <- !is.na(surv.iso) & surv.iso >= surv.cut
    upper <- if (any(upper.idx)) max(fit.times[upper.idx]) else max(fit.times)

    out <- c(lower, upper)
    out[1] <- max(out[1], min(fit.times))
    out[2] <- min(out[2], max(fit.times))
    if (out[1] >= out[2]) out <- c(min(fit.times), max(fit.times))
    return(out)
}

.np_contrast_band_endpts <- function(event.times, surv.0.iso, surv.1.iso,
                                     fit.times, uniform.cutpoint) {
    event.times <- event.times[is.finite(event.times)]
    if (length(event.times) == 0) return(c(min(fit.times), max(fit.times)))

    lower <- as.numeric(stats::quantile(event.times, uniform.cutpoint[1], na.rm = TRUE))
    surv.cut <- 1 - uniform.cutpoint[2]
    upper.idx <- (!is.na(surv.0.iso) & surv.0.iso >= surv.cut) |
        (!is.na(surv.1.iso) & surv.1.iso >= surv.cut)
    upper <- if (any(upper.idx)) max(fit.times[upper.idx]) else max(fit.times)

    out <- c(lower, upper)
    out[1] <- max(out[1], min(fit.times))
    out[2] <- min(out[2], max(fit.times))
    if (out[1] >= out[2]) out <- c(min(fit.times), max(fit.times))
    return(out)
}

.np_trim_uniform_window <- function(uniform.window, fit.times) {
    out <- c(max(uniform.window[1], min(fit.times)),
             min(uniform.window[2], max(fit.times)))
    if (!any(fit.times >= out[1] & fit.times <= out[2])) {
        stop("No `fit.times` fall inside `uniform.window`.")
    }
    return(out)
}

.np_uniform_test <- function(result, time, event, uniform.cutpoint=c(0.01, 0.99),
                             uniform.window=NULL, conf.level=.95, theta=0) {
    surv.0 <- .np_get_surv_object(result, trt = 0, isotonize = TRUE)
    surv.1 <- .np_get_surv_object(result, trt = 1, isotonize = TRUE)
    if (is.null(uniform.window)) {
        band.end.pts <- .np_contrast_band_endpts(time[event == 1],
                                                 surv.0$surv.iso, surv.1$surv.iso,
                                                 result$fit.times, uniform.cutpoint)
    } else {
        band.end.pts <- .np_trim_uniform_window(uniform.window, result$fit.times)
    }
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
