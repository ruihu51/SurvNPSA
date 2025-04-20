# estimate the identifiable components from the observed distribution
# input: observed data, nuisance estimators
# output: estimators and IF

.get.obs.comps <- function(time, event, treat, result, psi.type, verbose=TRUE){

    fit.times <- result$fit.times

    # estimate theta_n^obs = surv.1 - surv.0
    if(verbose) message("Estimating E[S(t|0,W)]...")
    surv.0 <- .get.survival(Y=time, Delta=event, A=1-treat, fit.times=fit.times,
                            eval.times=result$nuisance$eval.times,
                            S.hats=result$nuisance$event.pred.0,
                            G.hats=result$nuisance$cens.pred.0,
                            g.hats=1-result$nuisance$prop.pred)
    result$IF.vals.0 <- surv.0$IF.vals

    if(verbose) message("Estimating E[S(t|1,W)]...")
    surv.1 <- .get.survival(Y=time, Delta=event, A=treat, fit.times=fit.times,
                            eval.times=result$nuisance$eval.times,
                            S.hats=result$nuisance$event.pred.1,
                            G.hats=result$nuisance$cens.pred.1,
                            g.hats=result$nuisance$prop.pred)
    result$IF.vals.1 <- surv.1$IF.vals

    result$obs.comps.df <- data.frame(time=fit.times,
                                      theta.obs=surv.1$surv - surv.0$surv,
                                      surv.1=surv.1$surv,
                                      surv.0=surv.0$surv)
    result$IF.vals.theta.obs <- result$IF.vals.1 - result$IF.vals.0

    # estimate psi_n(t) = E{1(T>t)-S(t|A,W)}^2
    if(verbose) message("Estimating psi...")
    psi.rst <- .estimate.psi(Y=time, Delta=event, fit.times=fit.times,
                             eval.times=result$nuisance$eval.times,
                             S.hats=result$nuisance$event.pred,
                             G.hats=result$nuisance$cens.pred,
                             psi.type)
    result$obs.comps.df$psi <- psi.rst$psi.est
    result$obs.comps.df$gamma <- psi.rst$gamma.est
    result$IF.vals.psi <- psi.rst$IF.vals.psi

    # estimate tau_n = E(alpha_s(A,W)^2)
    if(verbose) message("Estimating tau...")
    tau.rst <- .estimate.tau(A=treat, g.hats=result$nuisance$prop.pred)
    result$tau <- tau.rst$est
    result$IF.vals.tau <- tau.rst$IF.vals.tau

    return(result)
}

# .get.obs.rmst <- function(t, result){
#     theta.obs <- result$obs.comps.df$theta.obs
#     theta.obs.func <- function(t0){
#         k <- min(which(fit.times >= t0))
#         return(theta.obs[k])
#     }
#
#     result$rmst.obs.est <- integrate(Vectorize(theta.obs.func), lower=0, upper=t)$value
#
#     .get.rmst.IF <- function(IF.vals.theta.obs, t){
#         theta.obs.IF.func <- function(t0){
#             k <- min(which(fit.times >= t0))
#             return(IF.vals.theta.obs[k])
#         }
#
#         IF.vals.rmst.obs <- tryCatch({
#             integrate(Vectorize(theta.obs.IF.func), lower=0, upper=t, subdivisions = 1000L)$value
#         }, error = function(e) {
#             message("An error occurred: ", e$message)
#             NA
#         })
#
#         return(IF.vals.rmst.obs)
#
#     }
#
#     result$IF.vals.rmst.obs <- apply(result$IF.vals.theta.obs, 1, .get.rmst.IF, t=t)
#     return(result)
# }


.get.survival <- function(Y, Delta, A, fit.times, eval.times, S.hats, G.hats, g.hats, isotonize=TRUE) {
    fit.times <- fit.times[fit.times > 0]
    n <- length(Y)
    ord <- order(eval.times)
    eval.times <- eval.times[ord]
    S.hats <- S.hats[,ord]
    G.hats <- G.hats[,ord]

    int.vals <- t(sapply(1:n, function(i) {
        vals <- diff(1/S.hats[i,])* 1/ G.hats[i,-ncol(G.hats)]
        if(any(eval.times[-1] > Y[i])) vals[eval.times[-1] > Y[i]] <- 0
        c(0,cumsum(vals))
    }))
    S.hats.Y <- sapply(1:n, function(i) stepfun(eval.times, c(1,S.hats[i,]), right = FALSE)(Y[i]))
    G.hats.Y <- sapply(1:n, function(i) stepfun(eval.times, c(1,G.hats[i,]), right = TRUE)(Y[i]))
    IF.vals <- matrix(NA, nrow=n, ncol=length(fit.times))
    surv <- rep(NA, length(fit.times))
    for(t0 in fit.times) {
        k <- min(which(eval.times >= t0))
        S.hats.t0 <- S.hats[,k]
        inner.func.1 <- ifelse(Y <= t0 & Delta == 1, 1/(S.hats.Y * G.hats.Y), 0 )
        inner.func.2 <- int.vals[,k]
        if.func <- as.numeric(A == 1) * S.hats.t0 * ( -inner.func.1 + inner.func.2) / g.hats + S.hats.t0
        k1 <- which(fit.times == t0)
        surv[k1] <- mean(if.func)
        IF.vals[,k1] <- if.func - surv[k1]
    }
    res <- list(times=fit.times, surv=pmin(1,pmax(0,surv)), IF.vals=IF.vals)
    if(isotonize) {
        res$surv.iso <- NA
        res$surv.iso[!is.na(res$surv)] <- 1 - isoreg(res$times[!is.na(res$surv)], 1-res$surv[!is.na(res$surv)])$yf
    }

    return(res)
}

.estimate.psi <- function(Y, Delta, fit.times, eval.times, S.hats, G.hats, psi.type){
    fit.times <- fit.times[fit.times > 0] # t0 strictly larger than 0
    n <- length(Y)
    ord <- order(eval.times)
    eval.times <- eval.times[ord]
    S.hats <- S.hats[,ord]
    G.hats <- G.hats[,ord]

    int.vals <- t(sapply(1:n, function(i) {
        vals <- diff(1/S.hats[i,])* 1/ G.hats[i,-ncol(G.hats)]
        if(any(eval.times[-1] > Y[i])) vals[eval.times[-1] > Y[i]] <- 0
        c(0,cumsum(vals))
    }))

    S.hats.Y <- sapply(1:n, function(i) stepfun(eval.times, c(1,S.hats[i,]), right = FALSE)(Y[i]))
    G.hats.Y <- sapply(1:n, function(i) stepfun(eval.times, c(1,G.hats[i,]), right = TRUE)(Y[i]))

    IF.vals.psi <- matrix(NA, nrow=n, ncol=length(fit.times))
    IF.vals.gamma <- matrix(NA, nrow=n, ncol=length(fit.times))
    psi.est <- rep(NA, length(fit.times))
    gamma.est <- rep(NA, length(fit.times))

    for(t0 in fit.times) {
        k <- min(which(eval.times >= t0))
        S.hats.t0 <- S.hats[,k]
        inner.func.1 <- ifelse(Y <= t0 & Delta == 1, 1/(S.hats.Y * G.hats.Y), 0 )
        inner.func.2 <- int.vals[,k]
        if.func.psi <- (1 - 2*S.hats.t0) * S.hats.t0 * ( -inner.func.1 + inner.func.2) +
            S.hats.t0 * (1 - S.hats.t0)
        k1 <- which(fit.times == t0)
        psi.est[k1] <- mean(if.func.psi)
        IF.vals.psi[,k1] <- if.func.psi - psi.est[k1]
        # plug-in estimator
        if (psi.type=="plug.in") {
            psi.est.plug.in <- mean(pmax(0,S.hats.t0)*(1-pmax(0,S.hats.t0)))
            psi.est[k1] <- psi.est.plug.in
        } else if (psi.type=="hybrid") {
            psi.est.plug.in <- mean(pmax(0,S.hats.t0)*(1-pmax(0,S.hats.t0)))
            psi.est[k1] <- ifelse(psi.est[k1]>0, psi.est[k1],
                                  psi.est.plug.in)
        } else {
            psi.est[k1] <- mean(if.func.psi)
        }

        if.func.gamma <- S.hats.t0^2 * 2 * ( -inner.func.1 + inner.func.2) + S.hats.t0^2
        gamma.est[k1] <- mean(if.func.gamma)
        IF.vals.gamma[,k1] <- if.func.gamma - gamma.est[k1]
    }
    res <- list(times=fit.times, psi.est=psi.est, IF.vals.psi=IF.vals.psi,
                gamma.est=gamma.est, IF.vals.gamma=IF.vals.gamma)

    return(res)
}


.estimate.tau <- function(A, g.hats){
    # estimate tau=E[\alpha_s^2]
    if.func.tau <- 2/(g.hats*(1-g.hats)) - ((A-g.hats)/(g.hats*(1-g.hats)))^2
    tau.est <- mean(if.func.tau)
    IF.vals.tau <- if.func.tau - tau.est
    res <- list(est=tau.est, IF.vals.tau=IF.vals.tau)

    return(res)

}



