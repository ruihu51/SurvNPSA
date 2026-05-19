.get.nuisances.est <- function(time, event, treat, confounders,
                               fit.times=sort(unique(time[time > 0 & time < max(time[event == 1])])),
                               fit.treat=c(0,1), nuisance.options = list(),
                               verbose=TRUE){
    .args <- mget(names(formals()), sys.frame(sys.nframe()))

    nuis <- do.call("CFsurvival.nuisance.options", nuisance.options)

    n <- length(time)

    if(sum(1-event) == 0) {
        message("No censored events; unanticipated errors may occur.")
        nuis$cens.pred.0 <- nuis$cens.pred.1 <- matrix(1, nrow=n, ncol=k)
    }
    if(sum(event) == 0) {
        stop("No uncensored events; cannot perform estimation.")
    }

    if(!nuis$cross.fit) {
        if(!is.null(nuis$folds)) {
            message("nuisance.options$cross.fit set to FALSE, but V > 1 or folds specified. Cross-fitting will not be performed.")
        }
        nuis$V <- 1
        nuis$folds <- rep(1, n)
    }
    if(!is.null(nuis$folds)) {
        nuis$folds <- as.numeric(factor(nuis$folds))
        nuis$V <- length(unique(nuis$folds))
    }
    if(nuis$cross.fit & is.null(nuis$folds)) {
        event.0 <- which(event == 0)
        event.1 <- which(event == 1)
        folds.0 <- sample(rep(1:nuis$V, length = length(event.0)))
        folds.1 <- sample(rep(1:nuis$V, length = length(event.1)))
        nuis$folds <- rep(NA, n)
        nuis$folds[event.0] <- folds.0
        nuis$folds[event.1] <- folds.1
    }

    if(is.null(nuis$verbose)) nuis$verbose <- FALSE

    .check.input(time=time, event=event, treat=treat, confounders=confounders, fit.times=fit.times, fit.treat=fit.treat, nuisance.options=nuis, conf.band=conf.band, conf.level=conf.level, contrasts=contrasts, verbose=verbose)

    if(any(fit.times < 0)) {
        fit.times <- fit.times[fit.times > 0]
        message("fit.times < 0 removed.")
    }
    if(any(fit.times == 0)) fit.times <- fit.times[fit.times > 0]
    if(any(fit.times > max(time[event == 1]))) {
        fit.times <- fit.times[fit.times <= max(time[event == 1])]
        message("fit.times > max(time[event == 1]) removed.")
    }

    if(is.null(nuis$eval.times)) nuis$eval.times <- sort(unique(c(0,time[time > 0 & time <= max(fit.times)], max(fit.times))))
    # if(is.null(nuis$eval.times)) nuis$eval.times <- fit.times
    k <- length(nuis$eval.times)

    confounders <- as.data.frame(confounders)
    result <- list(fit.times=fit.times)

    #### ESTIMATE PROPENSITY ####

    if(is.null(nuis$prop.pred)) {
        if(verbose) message("Estimating propensities...")
        nuis$prop.pred <- rep(NA, n)
        if(nuis$V > 1) {
            if(nuis$save.nuis.fits) result$prop.fits <- vector(mode='list',length=nuis$V)
            for(v in 1:nuis$V) {
                if(verbose) message(paste("Fold ", v, "..."))
                train <- nuis$fold != v
                test <- nuis$fold == v
                prop.fit <- .estimate.propensity(A=treat[train], W=confounders[train,,drop=FALSE], newW=confounders[test,, drop=FALSE], SL.library=nuis$prop.SL.library, fit.treat=fit.treat, prop.trunc=nuis$prop.trunc, save.fit = nuis$save.nuis.fits, verbose = nuis$verbose)
                nuis$prop.pred[test] <- prop.fit$prop.pred
                if(nuis$save.nuis.fits) result$prop.fits[[v]] <- prop.fit$prop.fit
            }
        } else {
            prop.fit <- .estimate.propensity(A=treat, W=confounders, newW=confounders, SL.library=nuis$prop.SL.library, fit.treat=fit.treat, prop.trunc=nuis$prop.trunc,  save.fit = nuis$save.nuis.fits, verbose = FALSE)
            nuis$prop.pred <- prop.fit$prop.pred
            if(nuis$save.nuis.fits) {
                result$prop.fit <- prop.fit$prop.fit
            }
        }
    }

    #### ESTIMATE CONDITIONAL SURVIVALS ####

    if((1 %in% fit.treat & (is.null(nuis$event.pred.1) | is.null(nuis$cens.pred.1))) | (0 %in% fit.treat & (is.null(nuis$event.pred.0) | is.null(nuis$cens.pred.0)))) {
        if(verbose) message("Estimating conditional survivals...")
        if(0 %in% fit.treat & is.null(nuis$event.pred.0)) {
            do.event.pred.0 <- TRUE
            nuis$event.pred.0 <- matrix(NA, nrow=n, ncol=k)
        } else do.event.pred.0 <- FALSE
        if(1 %in% fit.treat & is.null(nuis$event.pred.1)) {
            do.event.pred.1 <- TRUE
            nuis$event.pred.1 <- matrix(NA, nrow=n, ncol=k)
        } else do.event.pred.1 <- FALSE
        if(0 %in% fit.treat & is.null(nuis$cens.pred.0)) {
            do.cens.pred.0 <- TRUE
            nuis$cens.pred.0 <- matrix(NA, nrow=n, ncol=k)
        } else do.cens.pred.0 <- FALSE
        if(1 %in% fit.treat & is.null(nuis$cens.pred.1)) {
            do.cens.pred.1 <- TRUE
            nuis$cens.pred.1 <- matrix(NA, nrow=n, ncol=k)
        } else do.cens.pred.1 <- FALSE

        # S(t|A,W) and G(t|A,W)
        nuis$event.pred <- matrix(NA, nrow=n, ncol=k)
        nuis$cens.pred <- matrix(NA, nrow=n, ncol=k)

        if(nuis$V > 1) {
            if(nuis$save.nuis.fits) result$surv.fits <- vector(mode='list',length=nuis$V)
            nuis$event.coef <- nuis$cens.coef <- list()
            for(v in 1:nuis$V) {
                if(verbose) message(paste("Fold ", v, "..."))
                train <- nuis$fold != v
                test <- nuis$fold == v
                surv.fit <- .estimate.conditional.survival(Y=time[train], Delta=event[train], A=treat[train], newA=treat[test],
                                                           W=confounders[train,, drop=FALSE], newW=confounders[test,, drop=FALSE],
                                                           event.SL.library=nuis$event.SL.library,
                                                           fit.times=nuis$eval.times, fit.treat=fit.treat,
                                                           cens.SL.library=nuis$cens.SL.library, survSL.control=nuis$survSL.control,
                                                           survSL.cvControl = nuis$survSL.cvControl, cens.trunc=nuis$cens.trunc,
                                                           save.fit = nuis$save.nuis.fits, verbose = nuis$verbose)
                if(do.event.pred.0) nuis$event.pred.0[test,] <- surv.fit$event.pred.0
                if(do.event.pred.1) nuis$event.pred.1[test,] <- surv.fit$event.pred.1
                if(do.cens.pred.0) nuis$cens.pred.0[test,] <- surv.fit$cens.pred.0
                if(do.cens.pred.1) nuis$cens.pred.1[test,] <- surv.fit$cens.pred.1

                # S(t|A,W) and G(t|A,W)
                nuis$event.pred[test,] <- surv.fit$event.pred
                nuis$cens.pred[test,] <- surv.fit$cens.pred

                if(nuis$save.nuis.fits) result$surv.fits[[v]] <- surv.fit$surv.fit
                nuis$event.coef[[v]] <- surv.fit$event.coef
                nuis$cens.coef[[v]] <- surv.fit$cens.coef
            }
        } else {
            surv.fit <- .estimate.conditional.survival(Y=time, Delta=event, A=treat, W=confounders, newW=confounders, event.SL.library=nuis$event.SL.library, fit.times=nuis$eval.times, fit.treat=fit.treat, cens.SL.library=nuis$cens.SL.library, survSL.control=nuis$survSL.control, survSL.cvControl = nuis$survSL.cvControl, cens.trunc=nuis$cens.trunc, save.fit = nuis$save.nuis.fits,  verbose = nuis$verbose)
            if(do.event.pred.0) nuis$event.pred.0 <- surv.fit$event.pred.0
            if(do.event.pred.1) nuis$event.pred.1 <- surv.fit$event.pred.1
            if(do.cens.pred.0) nuis$cens.pred.0 <- surv.fit$cens.pred.0
            if(do.cens.pred.1) nuis$cens.pred.1 <- surv.fit$cens.pred.1

            # S(t|A,W) and G(t|A,W)
            nuis$event.pred <- surv.fit$event.pred
            nuis$cens.pred <- surv.fit$cens.pred

            if(nuis$save.nuis.fits) result$surv.fit <- surv.fit$surv.fit
            nuis$event.coef <- surv.fit$event.coef
            nuis$cens.coef <- surv.fit$cens.coef
        }
    }

    result$nuisance <- nuis

    return(result)
}


CFsurvival.nuisance.options <- function(cross.fit = TRUE, V = ifelse(cross.fit, 10, 1), folds = NULL, eval.times = NULL,
                                        event.SL.library = lapply(c("survSL.km", "survSL.coxph", "survSL.expreg", "survSL.weibreg", "survSL.loglogreg", "survSL.rfsrc"), function(alg) c(alg, "All")),
                                        event.pred.0 = NULL, event.pred.1 = NULL, event.pred = NULL,
                                        cens.SL.library = lapply(c("survSL.km", "survSL.coxph", "survSL.expreg", "survSL.weibreg", "survSL.loglogreg", "survSL.rfsrc"), function(alg) c(alg, "All")),
                                        cens.trunc=0, cens.pred.0 = NULL, cens.pred.1 = NULL, cens.pred = NULL,
                                        survSL.control = list(initWeightAlg = "survSL.rfsrc", verbose=FALSE), survSL.cvControl = list(V = 10), save.nuis.fits = FALSE,
                                        prop.SL.library = lapply(c("SL.mean", "SL.glm", "SL.earth", "SL.ranger", "SL.xgboost"), function(alg) c(alg, "All")), prop.trunc=0, prop.pred = NULL,
                                        verbose=FALSE) {
    list(cross.fit = cross.fit, V = V, folds = folds, eval.times = eval.times,
         event.SL.library = event.SL.library, event.pred.0 = event.pred.0, event.pred.1 = event.pred.1, event.pred = event.pred,
         cens.SL.library = cens.SL.library, cens.trunc=cens.trunc, cens.pred.0 = cens.pred.0, cens.pred.1 = cens.pred.1, cens.pred = cens.pred,
         survSL.control = survSL.control,
         survSL.cvControl = survSL.cvControl,
         save.nuis.fits = save.nuis.fits,
         prop.SL.library = prop.SL.library,  prop.trunc=prop.trunc, prop.pred = prop.pred, verbose=verbose)
}

.check.input <- function(time, event, treat, confounders, fit.times, fit.treat, nuisance.options, verbose, ...) {
    if(any(time < 0)) stop("Only non-negative event/censoring times allowed!")
    if(any(time == 0 & event == 1)) stop("Events at time zero not allowed.")
    if(any(!(event %in% c(0,1)))) stop("Event must be binary.")
    if(any(!(treat %in% c(0,1)))) stop("Treatment must be binary.")
    if(length(time) != length(event) | length(time) != length(treat)) stop("time, event, and treat must be n x 1 vectors")
    if(any(!(fit.treat %in% c(0,1)))) stop("fit.treat must be a subset of c(0,1).")
    if(!is.null(nuisance.options$event.pred.1)) {
        if(!all.equal(dim(nuisance.options$event.pred.1), c(length(time), length(nuisance.options$eval.times)))) {
            stop("event.pred must be an n x k matrix (n is number of observations, k is length of eval.times)")
        }
    }
    if(!is.null(nuisance.options$event.pred.0)) {
        if(!all.equal(dim(nuisance.options$event.pred.0), c(length(time), length(nuisance.options$eval.times)))) {
            stop("event.pred must be an n x k matrix (n is number of observations, k is length of eval.times)")
        }
    }
    if(!is.null(nuisance.options$event.pred)) {
        if(!all.equal(dim(nuisance.options$event.pred), c(length(time), length(nuisance.options$eval.times)))) {
            stop("event.pred must be an n x k matrix (n is number of observations, k is length of eval.times)")
        }
    }
    if(!is.null(nuisance.options$cens.pred.1)) {
        if(!all.equal(dim(nuisance.options$cens.pred.1), c(length(time), length(nuisance.options$eval.times)))) {
            stop("cens.pred must be an n x k matrix (n is number of observations, k is length of eval.times)")
        }
    }
    if(!is.null(nuisance.options$cens.pred.0)) {
        if(!all.equal(dim(nuisance.options$cens.pred.0), c(length(time), length(nuisance.options$eval.times)))) {
            stop("cens.pred must be an n x k matrix (n is number of observations, k is length of eval.times)")
        }
    }
    if(!is.null(nuisance.options$cens.pred)) {
        if(!all.equal(dim(nuisance.options$cens.pred), c(length(time), length(nuisance.options$eval.times)))) {
            stop("cens.pred must be an n x k matrix (n is number of observations, k is length of eval.times)")
        }
    }
    if(is.null(nuisance.options$event.SL.library)) {
        if(0 %in% fit.treat & is.null(nuisance.options$event.pred.0)) {
            stop("event.pred.0 must be provided if event.SL.library is not specified and 0 is a treatment of interest.")
        }
        if(1 %in% fit.treat & is.null(nuisance.options$event.pred.1)) {
            stop("event.pred.1 must be provided if event.SL.library is not specified and 1 is a treatment of interest.")
        }
    }
    if(is.null(nuisance.options$cens.SL.library)) {
        if(0 %in% fit.treat & is.null(nuisance.options$cens.pred.0)) {
            stop("cens.pred.0 must be provided if cens.SL.library is not specified and 0 is a treatment of interest.")
        }
        if(1 %in% fit.treat & is.null(nuisance.options$cens.pred.1)) {
            stop("cens.pred.1 must be provided if cens.SL.library is not specified and 1 is a treatment of interest.")
        }
    }
    # if(length(contrasts) > 0 & !identical(sort(fit.treat), c(0,1))) {
    #     warning("contrast specified but both treatment regimens not requested -- contrasts will not be provided. Re-run with fit.treat = c(0,1) for survival contrasts.")
    # }
    if(any(is.na(time) | is.na(event))) {
        stop("Missing time or event detected; missing data not allowed.")
    }
    if(any(is.na(confounders) | is.na(treat))) {
        warning("Missing confounders or exposures detected.")
    }
}

.estimate.propensity <- function(A, W, newW, SL.library, fit.treat, prop.trunc, save.fit, verbose) {
    ret <- list()
    library(SuperLearner)

    if ((length(SL.library)==1) & (SL.library[[1]][1]=="SL.gam.custom")){
        message("Using SL.gam simplified version for simulation")
        cat("start:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
        prop.fit <- SL.gam.custom(Y=A, X=W, newX=newW)
        ret$prop.pred <- c(prop.fit$pred)
        cat("end:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
    } else {
        prop.fit <- SuperLearner(Y=A, X=W, newX=newW, family='binomial',
                                 SL.library=SL.library, method = "method.NNLS", obsWeights = rep(1, length(A)), verbose = verbose)
        ret$prop.pred <- c(prop.fit$SL.predict)
    }

    if(save.fit) ret$prop.fit <- prop.fit
    # }
    if(1 %in% fit.treat) ret$prop.pred <- pmax(ret$prop.pred, prop.trunc)
    if(0 %in% fit.treat) ret$prop.pred <- pmin(ret$prop.pred, 1-prop.trunc)
    return(ret)
}

.estimate.conditional.survival <- function(Y, Delta, A, W, newA, newW, fit.times, fit.treat, event.SL.library, cens.SL.library, survSL.control, survSL.cvControl, cens.trunc, verbose, save.fit) {
    ret <- list(fit.times=fit.times)
    AW <- cbind(A, W)
    if(0 %in% fit.treat & 1 %in% fit.treat) {
        newAW <- rbind(cbind(A=0, newW), cbind(A=1, newW))
    } else {
        newAW <- cbind(A=fit.treat, newW)
    }

    # S(t|A,W) and G(t|A,W)
    newAW <- rbind(newAW, cbind(A=newA, newW))

    res <- require(survSuperLearner)
    if(!res) stop("Please install the package survSuperLearner via:\n devtools::install_github('tedwestling/survSuperLearner')")
    if(is.null(survSL.control)) survSL.control <- list(saveFitLibrary = save.fit)

    if ((length(event.SL.library)==1) & (event.SL.library[[1]][1]=="survSL.gam.custom")){
        message("Using survSL.gam simplified version for simulation")
        cat("start:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
        fit <- surv.gam(time = Y, event = Delta,  X = AW, newX = newAW, new.times = fit.times)
        cat("end:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
    } else {
        fit <- survSuperLearner(time = Y, event = Delta,  X = AW, newX = newAW, new.times = fit.times,
                                event.SL.library = event.SL.library, cens.SL.library = cens.SL.library,
                                verbose=verbose, control = survSL.control, cvControl = survSL.cvControl)
    }

    if(save.fit) ret$surv.fit <- fit

    pred_n <- nrow(newW)
    pred_total <- dim(fit$event.SL.predict)[1]

    if(pred_total!=pred_n*(length(fit.treat)+1)) {
        stop("The number of prediction results does not match the expected count based on the provided input parameters.")
    }

    if(0 %in% fit.treat) {
        ret$event.pred.0 <- fit$event.SL.predict[1:nrow(newW),]
        if(any(ret$event.pred.0 == 0)) ret$event.pred.0[ret$event.pred.0 == 0] <- min(ret$event.pred.0[ret$event.pred.0 > 0])
        ret$cens.pred.0 <- fit$cens.SL.predict[1:nrow(newW),]
        ret$cens.pred.0 <- pmax(ret$cens.pred.0, cens.trunc)
        if(any(ret$cens.pred.0 == 0)) ret$cens.pred.0[ret$cens.pred.0 == 0] <- min(ret$cens.pred.0[ret$cens.pred.0 > 0])
        if(1 %in% fit.treat) {
            ret$event.pred.1 <- fit$event.SL.predict[(pred_n+1):(pred_n*2),]
            if(any(ret$event.pred.1 == 0)) ret$event.pred.1[ret$event.pred.1 == 0] <- min(ret$event.pred.1[ret$event.pred.1 > 0])

            ret$cens.pred.1 <- fit$cens.SL.predict[(pred_n+1):(pred_n*2),]
            ret$cens.pred.1 <- pmax(ret$cens.pred.1, cens.trunc)
            if(any(ret$cens.pred.1 == 0)) ret$cens.pred.1[ret$cens.pred.1 == 0] <- min(ret$cens.pred.1[ret$cens.pred.1 > 0])
        }

        # S(t|A,W) and G(t|A,W)
        ret$event.pred <- fit$event.SL.predict[(pred_total-pred_n+1):pred_total,]
        if(any(ret$event.pred == 0)) ret$event.pred[ret$event.pred == 0] <- min(ret$event.pred[ret$event.pred > 0])
        ret$cens.pred <- fit$cens.SL.predict[(pred_total-pred_n+1):pred_total,]
        ret$cens.pred <- pmax(ret$cens.pred, cens.trunc)
        if(any(ret$cens.pred == 0)) ret$cens.pred[ret$cens.pred == 0] <- min(ret$cens.pred[ret$cens.pred > 0])

    } else {
        ret$event.pred.1 <- fit$event.SL.predict[1:nrow(newW),]
        if(any(ret$event.pred.1 == 0)) ret$event.pred.1[ret$event.pred.1 == 0] <- min(ret$event.pred.1[ret$event.pred.1 > 0])
        ret$cens.pred.1 <- fit$cens.SL.predict[1:nrow(newW),]
        ret$cens.pred.1 <- pmax(ret$cens.pred.1, cens.trunc)
        if(any(ret$cens.pred.1 == 0)) ret$cens.pred.1[ret$cens.pred.1 == 0] <- min(ret$cens.pred.1[ret$cens.pred.1 > 0])

        # S(t|A,W) and G(t|A,W)
        ret$event.pred <- fit$event.SL.predict[(pred_n+1):(pred_n*2),]
        if(any(ret$event.pred == 0)) ret$event.pred[ret$event.pred == 0] <- min(ret$event.pred[ret$event.pred > 0])
        ret$cens.pred <- fit$cens.SL.predict[(pred_n+1):(pred_n*2),]
        ret$cens.pred <- pmax(ret$cens.pred, cens.trunc)
        if(any(ret$cens.pred == 0)) ret$cens.pred[ret$cens.pred == 0] <- min(ret$cens.pred[ret$cens.pred > 0])
    }

    ret$event.coef <- fit$event.coef
    ret$cens.coef <- fit$cens.coef
    return(ret)
}

surv.gam <- function(time, event, X, newX, new.times, cts.num = 5, k=5){
    is_num   <- vapply(X, is.numeric, logical(1))
    nunique  <- vapply(X, function(v) length(unique(v)), integer(1))
    cont     <- names(X)[is_num & (nunique > cts.num)]
    disc     <- setdiff(colnames(X), cont)

    cont_terms <- if (length(cont)) paste0("s(", cont, ", k=", k, ")") else character(0)
    # cont_terms <- if (length(cont)) paste0("s(", cont, ")") else character(0)
    lin_terms  <- if (length(disc)) disc else character(0)
    rhs <- paste(c(cont_terms, lin_terms), collapse = " + ")
    if (rhs == "") rhs <- "1"

    gam.model <- as.formula(paste("time~", rhs))
    print(gam.model)
    fit.gam.event <- mgcv::gam(gam.model, family=mgcv::cox.ph(), data = X, method = "REML", weights=event, bs = 'ts', mgcv.tol = 1e-3)
    fit.gam.cens <- mgcv::gam(gam.model, family=mgcv::cox.ph(), data = X, method = "REML", weights=1-event, bs = 'ts', mgcv.tol = 1e-3)

    new.data <- data.frame(time=rep(new.times, each = nrow(newX)))
    for (col in names(newX)) new.data[[col]] <- rep(newX[[col]], length(new.times))

    event.predict <- predict(fit.gam.event, newdata=new.data, type="response", se=FALSE)
    event.predict <- matrix(event.predict, nrow = nrow(newX), ncol = length(new.times))

    cens.predict <- predict(fit.gam.cens, newdata=new.data, type="response", se=FALSE)
    cens.predict <- matrix(cens.predict, nrow = nrow(newX), ncol = length(new.times))

    out <- list(event.SL.predict = event.predict,
                cens.SL.predict = cens.predict)
    return(out)
}

SL.gam.custom <- function(Y, X, newX, cts.num = 5, k = 5, ...){
    is_num   <- vapply(X, is.numeric, logical(1))
    nunique  <- vapply(X, function(v) length(unique(v)), integer(1))
    cont     <- names(X)[is_num & (nunique > cts.num)]
    disc     <- setdiff(colnames(X), cont)

    cont_terms <- if (length(cont)) paste0("s(", cont, ", k=", k, ")") else character(0)
    lin_terms  <- if (length(disc)) disc else character(0)
    rhs <- paste(c(cont_terms, lin_terms), collapse = " + ")
    if (rhs == "") rhs <- "1"

    gam.model <- as.formula(paste("Y ~", rhs))
    print(gam.model)
    fit.gam.prop <- mgcv::gam(gam.model, family='binomial', data = cbind(Y, X))
    prop.predict <- predict(fit.gam.prop, newdata=newX, type="response", se=FALSE)
    out <- list(pred = prop.predict, fit = fit.gam.prop)
    class(out$fit) <- c("SL.gam.custom")
    return(out)
}
