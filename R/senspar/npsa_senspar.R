.simulate.senspar <- function(time, event, treat, confounders,
                              fit.times,
                              psi, tau, S.hat.obs,
                              num_drop = NULL, pct_drop = NULL, rep = 100, seed = 6741,
                              rmst = TRUE, fit.times.rmst = NULL, gamma = NULL, max_gap = NULL, tol=NULL){

  n_var <- ncol(confounders)

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

  if (length(fit.times) != length(psi)) {
    stop("The length of 'fit.times' must match the length of 'psi'.")
  }

  # if (rmst){
  #   gap <- max(diff(sort(fit.times[fit.times >= min(fit.times) & fit.times <= max(fit.times.rmst)])))
  #   if (gap > max_gap) {
  #     stop(sprintf("Error: fit.times is too sparse (gap = %.2f). Cannot guarantee RMST integration accuracy.", gap))
  #   }
  # }

  Gain.out.df <- data.frame(value = numeric(), t = numeric(), j = integer(), d = integer())
  Gain.out.phi.df <- data.frame(value = numeric(), t = numeric(), j = integer(), d = integer())
  Gain.trt.df <- data.frame(value = numeric(), j = integer(), d = integer())

  num_drop <- sort(unique(c(num_drop, 1, ceiling(0.5 * n_var))))
  num_drop <- num_drop[num_drop >= 1 & num_drop < n_var]

  for(d in num_drop){
    if (d==1){
        drop.index <- seq(1, n_var, by=1)
        J <- n_var
    } else {
        comb.drop <- combn(n_var, d)
        J <- ifelse(dim(comb.drop)[2] > rep, rep, dim(comb.drop)[2])
        set.seed(seed)
        drop.index <- comb.drop[,sample(1:dim(comb.drop)[2], J, replace = FALSE)]
    }


    for (j in 1:J){
      cat("d =", d, " j =", j, " Time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
      if (is.matrix(drop.index)){
        confounders.drop <- confounders[,-(drop.index[,j])]
      } else {
        confounders.drop <- confounders[,-(drop.index[j])]
      }

      result.sim.drop <- .get.nuisances.est(time = time,
                                            event = event,
                                            treat = treat,
                                            confounders = confounders.drop,
                                            fit.times = fit.times,
                                            nuisance.options = list(),
                                            verbose = FALSE)
      result.sim.drop <- .get.obs.comps(time=time, event=event, treat=treat,
                                        result=result.sim.drop,
                                        psi.type = "hybrid",
                                        verbose = FALSE)

      eval.idx <- sapply(result.sim.drop$fit.times, function(ft) {
        which.min(ifelse(result.sim.drop$nuisance$eval.times >= ft,
                         result.sim.drop$nuisance$eval.times, Inf))
      })

      fit.idx <- sapply(result.sim.drop$fit.times, function(ft) {
        which.min(ifelse(fit.times >= ft,
                         fit.times, Inf))
      })

      S.hat.obs.drop <- result.sim.drop$nuisance$event.pred

      V.g.matrix.psi <- colMeans((S.hat.obs - S.hat.obs.drop)^2)
      V.a.vector <- tau - result.sim.drop$tau

      Gain.out.matrix = pmax(0, V.g.matrix.psi[eval.idx] / psi[fit.idx]) # 1*t
      Gain.trt.vector = pmax(0, V.a.vector / tau) # 1*1

      cat(length(result.sim.drop$fit.times), "\n")
      Gain.out.df <- rbind(Gain.out.df, data.frame(C.Y.sq = Gain.out.matrix,
                                                   t = result.sim.drop$fit.times,
                                                   j = j,
                                                   d = d))

      Gain.trt.df <- rbind(Gain.trt.df, data.frame(C.A.sq = Gain.trt.vector / (1 - Gain.trt.vector),
                                                   j = j,
                                                   d = d))

      if (rmst){

        h.hat.obs <- t(apply(S.hat.obs, 1, function(row) {
          sapply(fit.times.rmst, function(t) .get.S.hat.int.vals(t, row, tol = tol))
        }))
        h.hat.obs.drop <- t(apply(S.hat.obs.drop, 1, function(row) {
          sapply(fit.times.rmst, function(t) .get.S.hat.int.vals(t, row, tol = tol))
        }))

        V.h.matrix.gamma <- colMeans((h.hat.obs - h.hat.obs.drop)^2)

        Gain.out.phi.matrix = pmax(0, V.h.matrix.gamma / gamma) # 1*t

        Gain.out.phi.df <- rbind(Gain.out.phi.df, data.frame(C.Y.sq = Gain.out.phi.matrix,
                                                             t = fit.times.rmst,
                                                             j = j,
                                                             d = d))
      }

    }
  }

  sens.df <- merge(Gain.out.df, Gain.trt.df, by = c("j", "d"))

  sens.df.mean <- sens.df %>%
    mutate(sens.par = C.Y.sq * C.A.sq) %>%
    group_by(d, t) %>%
    summarize(sens.par = mean(sens.par), .groups = "drop")

  senspar <- list(sens.df = sens.df,
                  sens.df.mean = sens.df.mean)

  if (rmst){
    sens.rmst.df <- merge(Gain.out.phi.df, Gain.trt.df, by = c("j", "d"))

    sens.rmst.df.mean <- sens.rmst.df %>%
      mutate(sens.par = C.Y.sq * C.A.sq) %>%
      group_by(d, t) %>%
      summarize(sens.par = mean(sens.par), .groups = "drop")

    senspar <- list(sens.df = sens.df,
                    sens.df.mean = sens.df.mean,
                    sens.rmst.df = sens.rmst.df,
                    sens.rmst.df.mean = sens.rmst.df.mean)
  }

  return(senspar)

}

.get.S.hat.int.vals <- function(t, S.hat.obs, tol){
  # step function
  theta.obs.func <- stepfun(c(0, result.sim.drop$nuisance$eval.times[-length(result.sim.drop$nuisance$eval.times)]),
                            c(1, S.hat.obs), right = FALSE)
  # integration
  hcubature(theta.obs.func, lowerLimit = c(0),
            upperLimit = c(t), tol=tol)$integral
}
