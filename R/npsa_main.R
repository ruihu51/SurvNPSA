
cat("start estimating nuisances:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
result <- .get.nuisances.est(time = time,
                             event = event,
                             treat = treat,
                             confounders = confounders,
                             fit.times = fit.times,
                             nuisance.options = list(),
                             verbose = TRUE)

cat("start estimating target:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
result <- .get.obs.comps(time=time, event=event, treat=treat, result=result,
                         psi.type = "hybrid",
                         verbose = TRUE)

cat("start estimating RMST:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
result <- .get.rmst.obs.comps(time = time,
                              event = event, 
                              result = result, 
                              fit.times.rmst = c(0.5, 0.7, 2), max_gap=0.2,
                              tol=0.01, tol1=0.01, tol2=0.01, 
                              gamma.type="hybrid", verbose=TRUE)

cat("start estimating RMST:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
senspar.df <- .simulate.senspar(time = time,
                                event = event,
                                treat = treat,
                                confounders = confounders,
                                fit.times = result$fit.times, 
                                psi = result$obs.comps.df$psi, 
                                tau = result$tau, 
                                S.hat.obs = result$nuisance$event.pred,
                                pct_drop = c(0.3, 0.7), 
                                rep = 10,
                                rmst = FALSE)

senspar.df <- .simulate.senspar(time = time,
                                event = event,
                                treat = treat,
                                confounders = confounders,
                                fit.times = result$fit.times, 
                                psi = result$obs.comps.df$psi, 
                                tau = result$tau, 
                                S.hat.obs = result$nuisance$event.pred,
                                pct_drop = c(0.3, 0.7), 
                                rep = 10,
                                rmst = TRUE,
                                fit.times.rmst = result$fit.times.rmst, 
                                gamma = result$gamma.est, 
                                tol = 0.01)

.report.bounds()

res.RV <- .report.RV(rv.time = c(0.3, 0.5, 1.2), 
                     result = result,
                     unif = FALSE)
summary(res.RV)

res.RV <- .report.RV(rv.time = c(0.3, 0.5, 1.2), 
                     result = result,
                     unif = TRUE)

.interpret.RV()


