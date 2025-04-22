npsa_surv <- function(time, event, treat, confounders, fit.times,
                      plot.times, rv.times,
                      rmst = TRUE,
                      sens.df.mean = NULL,
                      sens.rmst.df.mean = NULL,
                      var_names = NULL){
    if (!rmst){
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

        bounds.df <- .report.bounds(plot.times = c(0.5, 0.8, 1.2),
                                    result = result,
                                    rmst = FALSE)

        q.01 <- quantile(time[event == 1], .01)
        q.99 <- quantile(time[event == 1], .99)

        res.RV <- .report.RV(rv.time = c(0.3, 0.5, 1.2),
                             result = result,
                             unif = TRUE, q.01 = q.01, q.99 = q.99)

        summary(res.RV)

        if(is.null(sens.df.mean)){
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
        }

        bounds.df <- .report.bounds(plot.times = c(0.5, 0.8, 1.2),
                                    result = result,
                                    sens.df.mean = senspar.df$sens.df.mean,
                                    num_drop = NULL, pct_drop = c(0.3), n_var=10,
                                    rmst = FALSE)

        plot(bounds.df$bounds.df)

        interp.RV <- .interpret.RV(t0=1.2, res.RV,
                                   sens.df = senspar.df$sens.df,
                                   sens.df.mean = senspar.df$sens.df.mean,
                                   var_names, type = c("RV", "MIRV", "URV"))

        summary(interp.RV)

    } else{

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

        bounds.df <- .report.bounds(plot.times = c(0.5, 0.8, 1.2),
                                    result = result,
                                    rmst = rmst)

        if (!exists("sens.rmst.df.mean") || is.null(sens.rmst.df.mean)){
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
                                            max_gap = 0.2, tol=0.01)
        }

        bounds.df <- .report.bounds(plot.times = c(0.5, 0.8, 1.2),
                                    result = result,
                                    sens.df.mean = senspar.df$sens.df.mean,
                                    num_drop = NULL, pct_drop = c(0.3), n_var=10,
                                    rmst = rmst, sens.rmst.df.mean = senspar.df$sens.rmst.df.mean)

        plot(bounds.df$bounds.df)

    }
}







