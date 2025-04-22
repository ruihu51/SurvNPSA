#' Compute Effect Bounds (Internal)
#'
#' Internal utility function to compute the sensitivity analysis effect bounds.
#'
#' @param fit.times Numeric vector of times where bounds are evaluated.
#' @param theta.obs Numeric vector of observed differences.
#' @param psi Numeric vector of treatment assignment probabilities.
#' @param tau Truncation parameter.
#' @param sens.out Sensitivity parameter for outcome model.
#' @param sens.trt Sensitivity parameter for treatment model.
#' @param rho Correlation parameter.
#'
#' @return A list of lower and upper bounds at specified times.
#'
#' @keywords internal
.get.effect.bounds <- function(fit.times, theta.obs, psi, tau, sens.out, sens.trt, rho=1){
  # theta.obs <- c(1,2,3,5,7)
  # psi <- c(1,2,3,5,7)
  # sens.out <- c(1,2,3,5,7)
  # tau <- 0.1
  # sens.trt <- 0.3
  # fit.times <- c(1,2,3,5,6)
  # .get.bias.bounds(fit.times, theta.obs, psi, tau, sens.out, sens.trt)
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
                             rho=1, band.end.pts = c(0,Inf), conf.level=.95){
  # check input
  # times, bias.lower, bias.upper,
  # IF.vals.theta.obs, IF.vals.psi
  # sens.out length should match with each other
  fit.times <- effect.bounds$times
  effect.lower <- effect.bounds$effect.lower
  effect.upper <- effect.bounds$effect.upper
  sens.out <- effect.bounds$sens.out
  sens.trt <- effect.bounds$sens.trt

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

  se.lower <- sqrt(colMeans(IF.vals.effect.lower^2)) / sqrt(n) # 1*t
  se.upper <- sqrt(colMeans(IF.vals.effect.upper^2)) / sqrt(n) # 1*t

  # pointwise confidence intervals for effect bounds
  ptwise.bounds.lower <- effect.lower - qnorm(conf.level)*se.lower # z_{1-\alpha}
  ptwise.bounds.upper <- effect.upper + qnorm(conf.level)*se.upper

  # pointwise transformation
  trans.l <- log(1-pmin(effect.lower, 1))
  IF.trans.l <- IF.vals.effect.lower
  trans.u <- log(pmax(effect.upper, -1)+1)
  IF.trans.u <- IF.vals.effect.upper

  for(j in 1:length(effect.lower)) {
    IF.trans.l[,j] <- IF.vals.effect.lower[,j] * (1/(effect.lower[j]-1))
    IF.trans.u[,j] <- IF.vals.effect.upper[,j] * (1/(effect.upper[j]+1))
  }

  se.trans.l <- sqrt(colMeans(IF.trans.l^2))
  se.trans.u <- sqrt(colMeans(IF.trans.u^2))

  ptwise.trans.l <- 1 - exp(trans.l + qnorm(conf.level)*se.trans.l/sqrt(n)) # z_{1-\alpha}
  ptwise.trans.u <- exp(trans.u + qnorm(conf.level)*se.trans.u/sqrt(n)) - 1 # z_{1-\alpha}


  # uniform confidence band
  epsilon.l <- .estimate.limit.dist(IF.vals = IF.vals.effect.lower) # the output is boot*t
  # maxes <- replicate(1e4, max(rbind(rt(n, df = n - 1)/sqrt(n)) %*% IF.vals.effect.lower))
  # quantile.l <- quantile(maxes, conf.level)
  quantile.l <- quantile(apply(epsilon.l,1,max), conf.level)
  uniform.bounds.lower <- effect.lower - quantile.l / sqrt(n)

  epsilon.u <- .estimate.limit.dist(IF.vals = IF.vals.effect.upper)
  quantile.u <- quantile(apply(epsilon.u,1,min), 1-conf.level)
  uniform.bounds.upper <- effect.upper - quantile.u / sqrt(n)

  # uniform transformation

  # epsilon.trans.l <- .estimate.limit.dist(IF.vals = scale(IF.trans.l))
  # quantile.trans.l <- quantile(max(epsilon.trans.l), conf.level)
  # uniform.trans.l <- 1 - exp(trans.l + quantile.trans.l/sqrt(n))


  IF.trans.l.c <- IF.trans.l[,(fit.times >= band.end.pts[1] & fit.times <= band.end.pts[2])]
  # se.trans.l <- sqrt(colMeans(IF.trans.l^2))
  epsilon.trans.l <- .estimate.limit.dist(IF.vals = scale(IF.trans.l.c))
  quantile.trans.l <- quantile(apply(epsilon.trans.l,1,min), 1-conf.level)
  uniform.trans.l <- 1 - exp(trans.l - quantile.trans.l*se.trans.l/sqrt(n))
  uniform.trans.l[fit.times < band.end.pts[1] | fit.times > band.end.pts[2]] <- NA

  # epsilon.trans.u <- .estimate.limit.dist(IF.vals = scale(IF.trans.u))
  # quantile.trans.u <- quantile(max(epsilon.trans.u), conf.level)
  # uniform.trans.u <- exp(trans.u + quantile.trans.u/sqrt(n)) - 1

  IF.trans.u.c <- IF.trans.u[,(fit.times >= band.end.pts[1] & fit.times <= band.end.pts[2])]
  epsilon.trans.u <- .estimate.limit.dist(IF.vals = scale(IF.trans.u.c))
  quantile.trans.u <- quantile(apply(epsilon.trans.u,1,min), 1-conf.level)
  uniform.trans.u <- exp(trans.u - quantile.trans.u*se.trans.u/sqrt(n)) - 1
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

    bounds.df <- data.frame(
    times = bounds.conf.int$times,
    d = d
  )
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
    new_row <- as.data.frame(list(0, d, NA, 0, 0, 0, 0, 0, NA))
    names(new_row) <- names(bounds.df)

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
#' @param eval.times Numeric vector of times at which robustness values are evaluated.
#' @param rv Numeric vector of pointwise robustness values at each time.
#' @param q.01 Estimated lower bound (e.g., 1% survival quantile).
#' @param q.99 Estimated upper bound (e.g., 99% survival quantile).
#'
#' @return Numeric value representing the uniform robustness value.
#'
#' @keywords internal
.get.uniform.RV <- function(theta.obs, psi, tau,
                            IF.vals.theta.obs, IF.vals.psi, IF.vals.tau,
                            rho=1, conf.level=.95){

  # p-value under observed data
  # only proceed when p<0.05
  n <- nrow(IF.vals.theta.obs)
  test.stat <- n^(1/2)*sum(abs(theta.obs))
  epsilon <- .estimate.limit.dist(IF.vals = IF.vals.theta.obs)
  # dist.null <- apply(epsilon, 1, function(x) {max(abs(x))})
  dist.null <- rowSums(abs(epsilon)) # integration
  # dist.null <- replicate(1e4, sum(abs(rbind(rt(n, df = n - 1)/sqrt(n)) %*% IF.vals.theta.obs)))
  pvalue <- mean(dist.null > test.stat)
  cat("The p-value under no unobserved confounding is:", pvalue, "\n")

  if (pvalue < 1-conf.level) {

    message("Proceed to the test under unobserved confounding...")
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

      test.stat.sens <- max(n^(1/2)*max((theta.obs-sens.all*(sqrt(psi)*sqrt(tau)))),
                            n^(1/2)*max(-(theta.obs+sens.all*(sqrt(psi)*sqrt(tau)))))

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
      uniroot(pvalue.root, c(0.01,0.99), tol = 0.0001)$root
    }, error = function(e) {
      message("An error occurred: ", e$message)
      NA
    })
  } else {
    message("The null hypothesis that the observed effect is zero cannot be rejected. Sensitivity analysis will not proceed.")
    uniform.RV <- NA
  }

  return(uniform.RV)
}

.estimate.limit.dist.bound <- function(IF.vals.l, IF.vals.u, boot=10000){
  n <- nrow(IF.vals.l)
  t <- ncol(IF.vals.l)

  cov.matrix.l <- matrix(0, nrow = t, ncol = t)
  for (i in 1:n) {
    cov.matrix.l <- cov.matrix.l + outer(IF.vals.l[i,], IF.vals.l[i,])
  }
  cov.matrix.u <- matrix(0, nrow = t, ncol = t)
  for (i in 1:n) {
    cov.matrix.u <- cov.matrix.u + outer(IF.vals.u[i,], IF.vals.u[i,])
  }
  cov.matrix.l.u <- matrix(0, nrow = t, ncol = t)
  for (i in 1:n) {
    cov.matrix.l.u <- cov.matrix.l.u + outer(IF.vals.l[i,], IF.vals.u[i,])
  }
  cov.matrix.u.l <- t(cov.matrix.l.u)

  cov.matrix <- rbind(cbind(cov.matrix.l, cov.matrix.l.u), cbind(cov.matrix.u.l, cov.matrix.u))
  cov.matrix <- cov.matrix / n

  library(mvtnorm)
  epsilon <- rmvnorm(n=boot, mean=rep(0, 2*t), sigma = cov.matrix)
  epsilon.l <- epsilon[,1:t]
  epsilon.u <- epsilon[,(t+1):(2*t)]

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
#' @param eval.times Numeric vector of evaluation times.
#' @param theta.obs Numeric vector of observed treatment effect estimates at evaluation times.
#' @param effect.lower Numeric vector of lower bounds at evaluation times.
#' @param effect.upper Numeric vector of upper bounds at evaluation times.
#'
#' @return A numeric vector of pointwise robustness values corresponding to each evaluation time.
#'
#' @keywords internal
.get.RV <- function(t0, fit.times, theta.obs, psi, tau,
                    IF.vals.theta.obs, IF.vals.psi, IF.vals.tau,
                    rho=1, theta=0, conf.bounds=TRUE, transform=FALSE, conf.level=.95, verbose=TRUE){
  res.RV <- NULL

  k <- min(which(fit.times >= t0))
  theta.obs.t0 = theta.obs[k]
  psi.t0 = psi[k]

  v <- (theta.obs.t0 - theta)^2/(rho^2*psi.t0*tau)
  bounds.RV <- (- v + sqrt(v^2 + 4*v))/2

  res.RV <- list(t0=t0, bounds.RV=bounds.RV, theta=theta, bounds.int.RV=NA, conf.level=conf.level)

  if (conf.bounds){
    n <- dim(IF.vals.theta.obs)[1]

    # pointwise confidence intervals as function of sensitivity parameters
    bounds.senspar <- function(x, lower.b=TRUE){
      inner.func.1 <- (x/sqrt(1-x))*abs(rho)*(1/2)/(sqrt(psi.t0)*sqrt(tau))
      inner.func.2 <- (tau*IF.vals.psi[,k]+psi.t0*IF.vals.tau)*inner.func.1 # n*1
      if (lower.b){
        IF.vals.effect.lower <- IF.vals.theta.obs[,k] - inner.func.2
        se.lower <- sqrt(mean(IF.vals.effect.lower^2)) / sqrt(n)
        effect.lower.sp <- theta.obs.t0 - abs(rho)*sqrt(psi.t0)*sqrt(tau)*(x/sqrt(1-x))

        if (transform){
          trans.l <- log(1-pmin(effect.lower.sp, 1))
          IF.trans.l <- IF.vals.effect.lower * (1/(effect.lower.sp-1))
          se.trans.l <- sqrt(mean(IF.trans.l^2))
          return(1 - exp(trans.l + qnorm(conf.level)*se.trans.l/sqrt(n)) - theta)
        } else {
          return(effect.lower.sp - qnorm(conf.level)*se.lower - theta)
        }

      } else {
        IF.vals.effect.upper <- IF.vals.theta.obs[,k] + inner.func.2
        se.upper <- sqrt(mean(IF.vals.effect.upper^2)) / sqrt(n)
        effect.upper.sp <- theta.obs.t0 + abs(rho)*sqrt(psi.t0)*sqrt(tau)*(x/sqrt(1-x))
        if (transform){
          trans.u <- log(pmax(effect.upper.sp, -1)+1)
          IF.trans.u <- IF.vals.effect.upper * (1/(effect.upper.sp+1))
          se.trans.u <- sqrt(mean(IF.trans.u^2))
          return(exp(trans.u + qnorm(conf.level)*se.trans.u/sqrt(n)) - 1 - theta)
        } else {
          return(effect.upper.sp + qnorm(conf.level)*se.upper - theta)
        }

      }
    }

    # confidence intervals for theta.obs.t0 when sensitivity parameters are 0
    l.sp.0 <- bounds.senspar(0, lower.b = TRUE)
    u.sp.0 <- bounds.senspar(0, lower.b = FALSE)
    if (verbose) cat("pointwise confidence interval for theta.obs.t0 when sensitivity parameters are 0:", c(l.sp.0, u.sp.0), "\n")

    if (l.sp.0 <= 0 & u.sp.0 >= 0) {
      message("Pointwise CI cover the hypothesized value of theta; robustness values calculation for the lower/upper limit unnecessary.")
      bounds.int.RV <- 0
      lower.b <- NULL
    } else {
      lower.b <- ifelse(l.sp.0 > 0, TRUE, FALSE)
      bounds.int.RV <- tryCatch({
        uniroot(bounds.senspar, c(0,0.99), tol = 0.0001, lower.b = lower.b)$root
      }, error = function(e) {
        message("An error occurred: ", e$message)
        NA
      })
    }

    res.RV <- list(t0=t0, bounds.RV=bounds.RV, theta=theta, bounds.int.RV=bounds.int.RV, conf.level=conf.level, lower.b=lower.b)
  }

  return(res.RV)
}
