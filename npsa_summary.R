.report.RV <- function(rv.time, result, conf.level = .95, unif = TRUE) {
  res.list <- list()
  
  for (t0 in rv.time){
    res.RV <- .get.RV(
      t0,
      fit.times = result$fit.times,
      theta.obs = result$obs.comps.df$theta.obs,
      psi = result$obs.comps.df$psi,
      tau = result$tau,
      IF.vals.theta.obs = result$IF.vals.theta.obs,
      IF.vals.psi = result$IF.vals.psi,
      IF.vals.tau = result$IF.vals.tau,
      conf.level = 1 - (1 - conf.level)/2
    )
    
    res.list[[length(res.list) + 1]] <- list(
      t0 = res.RV$t0,
      theta = res.RV$theta,
      RV = res.RV$bounds.RV,
      MIRV = res.RV$bounds.int.RV,
      conf.level = res.RV$conf.level,
      lower.b = if (is.null(res.RV$lower.b)) NA else res.RV$lower.b
    )
  }
  
  res.table <- do.call(rbind, lapply(res.list, as.data.frame))
  rownames(res.table) <- NULL
  
  out <- list(res.table = res.table)
  
  if (unif){
    unif.RV <- .get.uniform.RV(
      theta.obs = result$obs.comps.df$theta.obs,
      psi = result$obs.comps.df$psi,
      tau = result$tau,
      IF.vals.theta.obs = result$IF.vals.theta.obs,
      IF.vals.psi = result$IF.vals.psi,
      IF.vals.tau = result$tau,
      conf.level = .95
    )
    out$unif.RV <- unif.RV
  }

  class(out) <- "reportRV"
  return(out)
}

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

.interpret.RV <- function(){
    (sp.point <- 0.01781872*0.01781872/(1-0.01781872))
    (sp.l.pw <- 0.01291943*0.01291943/(1-0.01291943))
    (sp.unif <- unif.RV*unif.RV/(1-unif.RV))
    
    col <- c("PERSONYRS", "respiratory_mort", "RF_PHYS_MODVIG_CURR2",
             "ENTRY_AGE", "SEX", "RACEI", "EDUCM", "HEALTH",
             "BMI_CUR1", "HEI2015_TOTAL_SCORE", "MPED_A_BEV_NOFOOD",
             "SMOKE_QUIT_DETAILED", "SMOKE_DOSE")
    
    sens.df %>%
      mutate(sens.par = C.Y.sq * C.A.sq) %>%
      filter(t==10 & d==1) %>%
      mutate(sig.point=sens.par>sp.point,
             sig.l.pw=sens.par>sp.l.pw,
             confounder=col[drop.index.list[[1]][j]+3])
    
    ### leave-d-out
    sens.df.mean %>%
      filter(t==10) %>%
      mutate(sig.point=sens.par>sp.point,
             sig.l.pw=sens.par>sp.l.pw)
    
    ### leave-half-out
    mean(sens.df %>%
           mutate(sens.par = C.Y.sq * C.A.sq) %>%
           filter(t==10 & d==5) %>%
           mutate(value = sens.par <= sp.point) %>%
           pull(value))
    
    mean(sens.df %>%
           mutate(sens.par = C.Y.sq * C.A.sq) %>%
           filter(t==10 & d==5) %>%
           mutate(value = sens.par <= sp.l.pw) %>%
           pull(value))
    
    # uniform
    ### leave-one-out
    sens.df %>%
      mutate(sens.par = C.Y.sq * C.A.sq) %>%
      group_by(d, j) %>%
      summarize(sens.par = max(sens.par)) %>%
      ungroup() %>%
      filter(d==1) %>%
      mutate(
        sig.unif=sens.par>sp.unif,
        confounder=col[drop.index.list[[1]][j]+3])
    
    ### leave-d-out
    
    sens.df %>%
      mutate(sens.par = C.Y.sq * C.A.sq) %>%
      group_by(d, j) %>%
      summarize(sens.par = max(sens.par)) %>%
      ungroup() %>%
      group_by(d) %>%
      summarize(sens.par = mean(sens.par)) %>%
      mutate(sig.unif=sens.par>sp.unif)
    
    
    ### leave-half-out
    mean(sens.df %>%
           mutate(sens.par = C.Y.sq * C.A.sq) %>%
           group_by(d, j) %>%
           summarize(sens.par = max(sens.par)) %>%
           ungroup() %>%
           filter(d==5) %>%
           mutate(value = sens.par <= sp.unif) %>%
           pull(value))
}

.report.bounds <- function(sens.df.mean, plot.time, result, num_drop = NULL, pct_drop = NULL, rmst=TRUE){
  
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
  
  bounds.df <- data.frame()
  
  for (d in num_drop){
    sens.out.true.input <- as.vector(sens.df.mean[sens.df.mean$d == d, "sens.par"])$sens.par
    sens.trt.true <- 1
    
    senspar.idx <- sapply(plot.time, function(x) {
      which(near(x, sens.df.mean$t[sens.df.mean$d == d]))
    })
    
    obs.est.idx <- sapply(plot.time, function(x) {
      which(near(x, result$fit.times))
    })
    
    # estimate lower and upper bound for true effect
    # theta.n.l and theta.n.u
    effect.bounds <-.get.effect.bounds(fit.times = result$fit.times[obs.est.idx],
                                       theta.obs = result$obs.comps.df$theta.obs[obs.est.idx],
                                       psi = result$obs.comps.df$psi[obs.est.idx],
                                       tau = result$tau,
                                       sens.out = sens.out.true.input[senspar.idx],
                                       sens.trt = sens.trt.true,
                                       rho = 1)
    
    # estimate pointwise confidence intervals and uniform confidence bands
    bounds.conf.int <- .bounds.confints(effect.bounds,
                                        psi = result$obs.comps.df$psi[obs.est.idx],
                                        tau = result$tau,
                                        IF.vals.theta.obs = result$IF.vals.theta.obs[,obs.est.idx],
                                        IF.vals.psi = result$IF.vals.psi[,obs.est.idx],
                                        IF.vals.tau = result$IF.vals.tau,
                                        conf.level=.975)
    
    df <- bounds2df(bounds.conf.int, theta.obs=result$obs.comps.df$theta.obs[obs.est.idx],
                           transform=TRUE, time.zero=TRUE)
    
    if (rmst){
      
    }
    
    bounds.df <- rbind(bounds.df, df)
  }
  
  
  
  
}
