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
                              band.end.pts=c(0,Inf), conf.level=.95) {
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
        if(any(!is.na(res$se))) {
            unif.vals <- .np_estimate.uniform.quantile(IF.vals[,!is.na(res$se), drop=FALSE],
                                                       conf.level, scale = FALSE)
            unif.quant <- unif.vals$quantile
            out$ew.sim.maxes <- unif.vals$maxes
            out$unif.ew.quant <- unif.quant
            res$unif.ew.lower <- pmax(est - unif.quant / sqrt(n), 0)
            res$unif.ew.upper <- pmin(est + unif.quant / sqrt(n), 1)

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
