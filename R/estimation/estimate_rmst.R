library(cubature)

.get.rmst.obs.comps <- function(time, event, result,
                                fit.times.rmst, max_gap=0.1,
                                tol=0.001, tol1=0.001, tol2=0.001,
                                gamma.type="hybrid", verbose=TRUE){

    if (any(fit.times.rmst > max(result$fit.times))) {
        message("Some fit.times.rmst > maximum observed event time - removed for RMST estimation.")
        fit.times.rmst <- fit.times.rmst[fit.times.rmst <= max(result$fit.times)]
    }

    eval.times.rmst <- result$fit.times

    gap <- max(diff(sort(eval.times.rmst[eval.times.rmst >= min(eval.times.rmst) & eval.times.rmst <= max(fit.times.rmst)])))
    if (gap > max_gap) {
        stop(sprintf(
            paste0(
                "Error: fit.times is too sparse (gap = %.2f).\n",
                "Cannot guarantee RMST integration accuracy.\n",
                "Please re-estimate nuisances using a denser fit.times grid with `.get.nuisances.est()`."
            ),
            gap
        ))
    }

    if(verbose) message("Estimating RMST differences...")
    result$rmst.obs <- sapply(fit.times.rmst, .get.obs.rmst.int.vals, result$obs.comps.df$theta.obs, tol = tol)

    result$IF.vals.rmst.obs <- t(apply(result$IF.vals.theta.obs, 1, function(row) {
        sapply(fit.times.rmst, function(t) .get.obs.rmst.int.vals(t, row, tol = tol))
    }))

    if(verbose) message("Estimating E[(min(T,t)-S(t|A,W))^2]...")
    gamma.rst <- .estimate.gamma(Y=time,
                                 Delta=event,
                                 fit.times=eval.times.rmst,
                                 eval.times=result$nuisance$eval.times,
                                 S.hats=result$nuisance$event.pred,
                                 G.hats=result$nuisance$cens.pred
    )

    if(verbose) message("Estimating E[(min(T,t)-S(t|A,W))^2] part 1...")
    gamma.p1 <- sapply(fit.times.rmst, .get.gamma.p1.int.vals, gamma.rst$gamma.p1.est, tol = tol1)
    gamma.p1.plug.in <- sapply(fit.times.rmst, .get.gamma.p1.int.vals, gamma.rst$gamma.p1.est.plug.in, tol = tol1)

    IF.vals.gamma.p1.int <- t(apply(gamma.rst$IF.vals.gamma.p1, 1, function(row) {
        sapply(fit.times.rmst, function(t) .get.gamma.p1.int.vals(t, row, tol = tol1))
    }))

    if(verbose) message("Estimating E[(min(T,t)-S(t|A,W))^2] part 2...")
    eval.times.area <- expand.grid(eval.times.rmst, eval.times.rmst)

    gamma.p2 <- sapply(fit.times.rmst, .get.gamma.p2.int.vals, gamma.rst$gamma.p2.est, tol = tol2, eval.times.area=eval.times.area)
    gamma.p2.plug.in <- sapply(fit.times.rmst, .get.gamma.p2.int.vals, gamma.rst$gamma.p2.est.plug.in, tol = tol2, eval.times.area=eval.times.area)

    IF.vals.gamma.p2.int <- t(apply(gamma.rst$IF.vals.gamma.p2, 1, function(row) {
        sapply(fit.times.rmst, function(t) .get.gamma.p2.int.vals(t, row, tol = tol2, eval.times.area=eval.times.area))
    }))

    gamma.est <- gamma.p1 - gamma.p2
    gamma.est.plug.in <- gamma.p1.plug.in - gamma.p2.plug.in

    if (gamma.type=="plug.in"){
        result$gamma.est <- gamma.est.plug.in
    } else if (gamma.type=="hybrid") {
        result$gamma.est <- ifelse(gamma.est>0, gamma.est, gamma.est.plug.in)
    } else {
        result$gamma.est <- gamma.est
    }

    result$IF.vals.gamma <- IF.vals.gamma.p1.int - IF.vals.gamma.p2.int

    result$fit.times.rmst <- fit.times.rmst

    return(result)
}


.get.obs.rmst.int.vals <- function(t, theta.obs, tol){
    # step function
    theta.obs.func <- stepfun(c(0, eval.times.rmst[-length(eval.times.rmst)]),
                              c(0, theta.obs), right = FALSE)
    # integration
    hcubature(theta.obs.func, lowerLimit = c(0),
              upperLimit = c(t), tol=tol)$integral
}

.estimate.gamma <- function(Y, Delta, fit.times, eval.times, S.hats, G.hats){
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

    # estimate gamma part 1 E[S(t|A,W)]
    # RETURN a function of fit.times
    IF.vals.gamma.p1 <- matrix(NA, nrow=n, ncol=length(fit.times))
    gamma.p1.est <- rep(NA, length(fit.times))
    gamma.p1.est.plug.in <- rep(NA, length(fit.times))

    for(t0 in fit.times) {
        k <- min(which(eval.times >= t0))
        S.hats.t0 <- S.hats[,k]
        inner.func.1 <- ifelse(Y <= t0 & Delta == 1, 1/(S.hats.Y * G.hats.Y), 0 )
        inner.func.2 <- int.vals[,k]
        if.func.gamma.p1 <- S.hats.t0 * (1 -inner.func.1 + inner.func.2)

        k1 <- which(fit.times == t0)
        gamma.p1.est.plug.in[k1] <- mean(pmax(0,S.hats.t0))
        gamma.p1.est[k1] <- mean(if.func.gamma.p1)
        IF.vals.gamma.p1[,k1] <- if.func.gamma.p1 - gamma.p1.est[k1]
    }

    # estimate gamma part 2 E[S(u|A,W)*S(v|A,W)]
    # RETURN a function of c(u,v)
    IF.vals.gamma.p2 <- matrix(NA, nrow=n, ncol=length(fit.times)^2)
    gamma.p2.est <- rep(NA, length(fit.times)^2)
    gamma.p2.est.plug.in <- rep(NA, length(fit.times)^2)

    fit.time.area <- expand.grid(fit.times, fit.times)
    for(row in 1:nrow(fit.time.area)) {
        u0 <- fit.time.area[row,1]
        v0 <- fit.time.area[row,2]
        u_k <- min(which(eval.times >= u0))
        v_k <- min(which(eval.times >= v0))
        S.hats.u0 <- S.hats[,u_k]
        S.hats.v0 <- S.hats[,v_k]
        inner.func.u1 <- ifelse(Y <= u0 & Delta == 1, 1/(S.hats.Y * G.hats.Y), 0 )
        inner.func.u2 <- int.vals[,u_k]
        inner.func.v1 <- ifelse(Y <= v0 & Delta == 1, 1/(S.hats.Y * G.hats.Y), 0 )
        inner.func.v2 <- int.vals[,v_k]
        if.func.gamma.double.p1 <- S.hats.v0 * S.hats.u0 * (-inner.func.u1 + inner.func.u2)
        if.func.gamma.double.p2 <- S.hats.u0 * S.hats.v0 * (-inner.func.v1 + inner.func.v2)
        if.func.gamma.double.p3 <- S.hats.u0 * S.hats.v0

        k1 <- which((fit.time.area[,1] == u0) & (fit.time.area[,2] == v0))
        gamma.p2.est.plug.in[k1] <- mean(pmax(0,S.hats.u0)*pmax(0,S.hats.v0))
        gamma.p2.est[k1] <- mean(if.func.gamma.double.p1 + if.func.gamma.double.p2 + if.func.gamma.double.p3)
        IF.vals.gamma.p2[,k1] <- if.func.gamma.double.p1 + if.func.gamma.double.p2 + if.func.gamma.double.p3 - gamma.p2.est[k1]
    }

    res_df <- list(gamma.p1.est=gamma.p1.est,
                   gamma.p1.est.plug.in=gamma.p1.est.plug.in,
                   IF.vals.gamma.p1=IF.vals.gamma.p1,
                   gamma.p2.est=gamma.p2.est,
                   gamma.p2.est.plug.in=gamma.p2.est.plug.in,
                   IF.vals.gamma.p2=IF.vals.gamma.p2)
}


.get.gamma.p1.int.vals <- function(t, gamma.p1.est, tol){
    # step function
    gamma.p1.func <- function(u){
        step_fun <- stepfun(c(0, eval.times.rmst[-length(eval.times.rmst)]),
                            c(1, gamma.p1.est), right = FALSE)
        return(step_fun(u)*2*u)
    }
    # integration
    hcubature(gamma.p1.func, lowerLimit = c(0),
              upperLimit = c(t), tol = tol)$integral
}


.get.gamma.p2.int.vals <- function(t, gamma.p2.est, tol, eval.times.area){
    gamma.p2.func <- function(u){
        # eval.times.area[min(which((eval.times.area[,1] >= 0.56) & (eval.times.area[,2] >= 0.78))),]
        return(gamma.p2.est[min(which((eval.times.area[,1] >= u[1]) & (eval.times.area[,2] >= u[2])))])

    }
    hcubature(gamma.p2.func, lowerLimit = c(0,0),
              upperLimit = c(t,t), tol = tol)$integral
}
