#*********************************************************************#
#*** JMDEM (Joint Mean and Dispersion Effects Model) Version 1.0.1 ***#
#*********************************************************************#
jmdem <- function(mformula, dformula, data, mfamily = gaussian, dfamily = Gamma, 
                  weights, subset, dev.type = c("deviance", "pearson"), 
                  moffset = NULL, doffset = NULL, mustart = NULL, phistart = NULL, 
                  betastart = NULL, lambdastart = NULL, hessian = TRUE, na.action, 
                  grad.func = TRUE, fit.method = "jmdem.fit", 
                  method = c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN", "Brent"), 
                  df.adj = FALSE, disp.adj = FALSE, full.loglik = FALSE, 
                  beta.first = TRUE, prefit = TRUE, mcontrasts = NULL, 
                  dcontrasts = NULL, control = list(...), 
                  minv.method = c("solve", "chol2inv", "ginv"), ...)
{
    ### Listing of all model arguments
    orgmatch <- match.call()
    arglist <- as.list(orgmatch)[-1]
    argnam <- names(arglist)
    
    ### Check whether some character inputs are correct
    dev.type <- match.arg(dev.type)
    method <- match.arg(method)
    minv.method <- match.arg(minv.method, several.ok = TRUE)
    
    ### Create the function argument lists of jmdem.prefit
    arg.jmdem.prefit <- as.list(args(jmdem.prefit))
    arg.jmdem.prefit <- arg.jmdem.prefit[!names(arg.jmdem.prefit) == "" & !names(arg.jmdem.prefit) == "..."]

    ### Create the function argument lists of jmdem.fit
    arg.jmdem.fit <- as.list(args(jmdem.fit))
    arg.jmdem.fit <- arg.jmdem.fit[!names(arg.jmdem.fit) == "" & !names(arg.jmdem.fit) == "..."]

    ### Create the function argument lists of glm.fit
    arg.glm.fit <- as.list(args(glm.fit))
    arg.glm.fit <- arg.glm.fit[!names(arg.glm.fit) == "" & !names(arg.glm.fit) == "..."]

    ### Create the function argument lists of loglik
    arg.loglik <- as.list(args(loglik))
    arg.loglik <- arg.loglik[!names(arg.loglik) == "" & !names(arg.loglik) == "..."]

    ### Create the function argument lists of optim
    arg.optim <- as.list(args(optim)) #Function optim
    arg.optim <- arg.optim[!names(arg.optim) == "" & !names(arg.optim) == "..."]

    ### Check distribution family
    mfamily <- check.family(family = mfamily)
    dfamily <- check.family(family = dfamily, type = "dispersion")

    ### Update the families in the extend call list with the adjusted format
    arglist$mfamily <- mfamily
    arglist$dfamily <- dfamily
    
    ### Update the model formula in the extend call list with the adjusted format
    arglist$mformula <- as.formula(mformula)
    arglist$dformula <- as.formula(dformula)

    ### Separate terms of both the mean and dispersion models
    mintercept <- attr(terms(arglist$mformula), "intercept") == 1
    dintercept <- attr(terms(arglist$dformula), "intercept") == 1
    mtermnam <- attr(terms(arglist$mformula), "term.labels")
    dtermnam <- attr(terms(arglist$dformula), "term.labels")

    ### Save the original dataset for security
    data.org <- data

    ### Dealing with missing values
    if (match("subset", names(arglist), 0L)) {data <- data[eval(arglist$subset, data),]}

    ### Dealing with missing values
    na.actionnam <- if (is.null(arglist$na.action)) {getOption("na.action")} else {deparse(substitute(na.action))}
    na.action <- get(na.actionnam)
    arglist$na.action <- na.action
    data <- do.call(na.actionnam, list(object = data))
    
    ### Extract temporary response variable
    temp.y <- data[[all.vars(arglist$mformula)[1]]]

    ### Update the data for the dispersion model if the model is nested
    if (length(which(dtermnam %in% "eta")))
    {
        data$eta <- if (length(which(argnam %in% "mustart"))) {mfamily$linkfun(mustart)} else {rep(mfamily$linkfun(mean(temp.y)), nrow(data))}
    }
    
    ### Update the data for the mean model if the model is nested
    if (length(which(mtermnam %in% "delta")))
    {
        data$delta <- if (length(which(argnam %in% "phistart"))) {dfamily$linkfun(phistart)} else {rep(var(temp.y) / mfamily$variance(mean(temp.y)), nrow(data))}
    }
    
    ### Update data frames for fitting
    arglist$data <- data

    ### Evaluate data frame
    mframe <- model.frame(arglist$mformula, data = data, na.action = na.action) #model frame for mean model
    dframe <- model.frame(arglist$dformula, data = data, na.action = na.action) #model frame for dispersion model
    mterms <- attr(mframe, "terms")
    dterms <- attr(dframe, "terms")
    if(identical(fit.method, "model.frame")) return(list(mean = mframe, dispersion = dframe))
    
    ### Create data frame
    y <- model.response(mframe) #response variable
    x <- model.matrix(eval(arglist$mformula), data = mframe, contrasts.arg = mcontrasts) #design matrix for mean model
    z <- model.matrix(eval(arglist$dformula), data = dframe, contrasts.arg = dcontrasts) #design matrix for dispersion model
    n <- nrow(data)
    
    ### Count number of parameters
    nmterms <- ncol(x)
    ndterms <- ncol(z)

    ### Data frame attribtues
    nmissing <- attr(data, "na.action")
    xlevels <- .getXlevels(attr(mframe, "terms"), mframe)
    zlevels <- .getXlevels(attr(dframe, "terms"), dframe)
    ynames <- names(y)
    xnames <- colnames(x)
    znames <- colnames(z)
    
    ### Response mean and dispersion
    mY <- mean(y)
    vY <- var(y)
    dY <- vY / mfamily$variance(mY)
    
    ### Check and extend the initial response mean in the call list
    if (length(which(argnam %in% "mustart"))) 
    {
        arglist$mustart <- rep(eval(arglist$mustart), length.out = n)
    } 
    else 
    {
        arglist$mustart <- rep(mY, length.out = n)
    }
    
    ### Check and extend the response dispersion in the call list
    if (length(which(argnam %in% "phistart"))) 
    {
        arglist$phistart <- rep(eval(arglist$phistart), length.out = n)
    } 
    else 
    {
        arglist$phistart <- rep(dY, length.out = n)
    }
    
    ### Check and extend the initial weights in the call list
    if (length(which(argnam %in% "weights"))) 
    {
        arglist$weights <- rep(eval(arglist$weights), length.out = n)
    } 
    else 
    {
        arglist$weights <- rep(1, length.out = n)
    }
    
    ### Check and extend call list by offset
    arglist$moffset <- moffset <- if (length(which(argnam %in% "moffset"))) {eval(arglist$moffset, envir = data)} else {rep(0L, length.out = n)} + if (is.null(model.offset(mframe))) {rep(0L, length.out = n)} else {model.offset(mframe)}
    arglist$doffset <- doffset <- if (length(which(argnam %in% "doffset"))) {eval(arglist$doffset, envir = data)} else {rep(0L, length.out = n)} + if (is.null(model.offset(dframe))) {rep(0L, length.out = n)} else {model.offset(dframe)}

    ### Set start values for beta and lambda if prefit is not preferred
    if (!is.null(prefit) & !prefit)
    {
        if (is.null(arglist$betastart)) {arglist$betastart <- if (mintercept) {setNames(c(mY, rep(0, length(xnames) - 1)), xnames)} else {setNames(rep(0, length(xnames)), xnames)}}
        if (is.null(arglist$lambdastart)) {arglist$lambdastart <- if (dintercept) {setNames(c(dY, rep(0, length(znames) - 1)), znames)} else {setNames(rep(0, length(znames)), znames)}}
    }

    ### Initiate dummy parameter vector for later checking
    npar <- nmterms + ndterms
    beta.index <- 1:nmterms
    lambda.index <- nmterms + (1:ndterms)
    
    ### Get system or user control parameters for optim
    if (!is.null(arglist$control))
    {
        cargs <- lwhich(eval(arglist$control), as.list(args(jmdem.control)))
    } 
    else 
    {
        cargs <- list()
    }
    arglist$control <- do.call("jmdem.control", cargs)
    null.approx <- arglist$control$null.approx
    if (length(arglist$control$parscale) < npar) {arglist$control$parscale <- rep(arglist$control$parscale, npar)}
    if (length(arglist$control$ndeps) < npar) {arglist$control$ndeps <- rep(arglist$control$ndeps, npar)}


    ### Get initial value of the model parameters
    prefarglist <- c(list(y = y, x = x, z = z, xnames = xnames, znames = znames, nmterms = nmterms, ndterms = ndterms, mintercept = mintercept, dintercept = dintercept), 
                     lwhich(arglist, arg.jmdem.prefit))
    if (!is.null(prefit) & !prefit)
    {
        if (!is.null(arglist$betastart) & !mintercept) {prefarglist$moffset <- rep(mY, n)}
        if (!is.null(arglist$lambdastart) & !dintercept) {prefarglist$doffset <- rep(dY, n)}
    }
    init.est <- do.call("jmdem.prefit", prefarglist)

    ### Create the list for final summary
    init.mean.fit <- list(coefficients = init.est$beta, linear.predictors = init.est$eta, fitted.values = init.est$mu)
    init.dispersion.fit <- list(coefficients = init.est$lambda, linear.predictors = init.est$delta, fitted.values = init.est$phi)

    ### Mean model initial estimates
    beta <- init.est$beta
    mu <- init.est$mu
    eta <- init.est$eta
    
    ### Dispersion model initial estimates
    lambda <- init.est$lambda
    phi <- init.est$phi
    delta <- init.est$delta
    
    ### Update the data frames for jmdem if necessary
    if (length(which(mtermnam %in% "delta"))) {x[, "delta"] <- delta}
    if (length(which(dtermnam %in% "eta"))) {z[, "eta"] <- eta}

    ### Model argument list for the jmdem
    marglist <- list(x = x, y = y, z = z, weights = arglist$weights, beta = beta, lambda = lambda, 
                     mu = mu, phi = phi, method = method, hessian = hessian, na.action = arglist$na.action, 
                     mintercept = mintercept, dintercept = dintercept)
    
    ### Function argument list for the jmdem.fit, loglik, optim
    farglist.jmdem.fit <- lwhich(arg.jmdem.fit, marglist, reverse = TRUE)
    farglist.loglik <- lwhich(arg.loglik, marglist, reverse = TRUE)
    farglist.optim <- lwhich(arg.optim, marglist, reverse = TRUE)

    ### Append all unique function arguments
    farglist <- c(farglist.jmdem.fit, farglist.loglik, farglist.optim)
    farglist <- farglist[!duplicated(names(farglist))]
    
    ### Append model arguments and function arguments
    farglist <- append(marglist, lwhich(arglist, farglist))

    ### Estimation
    est <- do.call("jmdem.fit", farglist)
    
    ### Finalise the data frames for summary
    eta <- est$mean.linear.predictors
    delta <- est$dispersion.linear.predictor
    if (length(which(mtermnam %in% "delta"))) 
    {
        mframe[, "delta"] <- delta
        x[, "delta"] <- delta
    }
    if (length(which(dtermnam %in% "eta"))) 
    {
        dframe[, "eta"] <- eta
        z[, "eta"] <- eta
    }

    ### Design matrix of the null mean model
    x0 <- if (mintercept) {x[, "(Intercept)", drop = FALSE]} else {as.matrix(rep(0, n))}

    ### Prefit null model
    xnames0 <- if (mintercept) {"(Intercept)"} else {as.character(0)}
    nmterms0 <- ncol(x0)
    betastart0 <- betastart[1]
    prefarglist0 <- c(list(y = y, x = x0, z = z, xnames = xnames0, znames = znames, nmterms = nmterms0, ndterms = ndterms, mintercept = mintercept, dintercept = dintercept), 
                      lwhich(arglist, arg.jmdem.prefit))
    prefarglist0$betastart <- betastart0
    init.mest0 <- do.call("jmdem.prefit", prefarglist0)

    ### Fit null model
    farglist0 <- farglist
    farglist0$x <- x0
    farglist0$beta <- init.mest0$beta
    farglist0$mu <- init.mest0$mu
    farglist0$lambda <- init.mest0$lambda
    farglist0$phi <- init.mest0$phi
    farglist0$control$parscale <- farglist0$control$parscale[c(1, lambda.index)]
    farglist0$control$ndeps <- farglist0$control$ndeps[c(1, lambda.index)]
    mest0 <- do.call("jmdem.fit", farglist0)

    ### Fit dispersion model for deviance and null deviance
    dfarglist <- list(x = z, y = est$residuals ^ 2, weights = rep(0.5, n), family = dfamily, intercept = dintercept, 
                      offset = doffset, control = list(maxit = farglist0$control$maxit, epsilon = farglist0$control$epsilon, trace = farglist0$control$fit.trace))
    dest <- do.call("glm.fit", dfarglist)

    ### Deviance of the null mean model
    null.deviance <- mest0$deviance
    df.null <- mest0$df.residual
    
    ### Deviance and null deviance of the dispersion model
    d.deviance <- dest$deviance
    d.df.residual <- dest$df.residual
    d.null.deviance <- dest$null.deviance
    d.df.null <- dest$df.null
    
    ### Coefficients of the null model
    beta.null <- mest0$beta
    lambda.null <- mest0$lambda
    
    ### Save user-defined matrix inverse method
    if (length(minv.method) > 1) {minv.method <- NULL}
    
    ### Summary list
    fit <- c(est, list(control = arglist$control, data = data.org, mean.model = mframe, 
                       dispersion.model = dframe, call = orgmatch, mean.formula = mformula, 
                       dispersion.formula = dformula, fit.method = fit.method, 
                       mean.offset = moffset, dispersion.offset = doffset, 
                       dispersion.deviance = d.deviance, dispersion.df.residual = d.df.residual,
                       null.deviance = null.deviance, df.null = df.null,
                       dispersion.null.deviance = d.null.deviance, dispersion.df.null = d.df.null,
                       beta.null = beta.null, lambda.null = lambda.null, 
                       dispersion.null = mest0$dispersion, residuals.null = mest0$residuals,
                       mustart = mustart, phistart = phistart, betastart = betastart, 
                       lambdastart = lambdastart, mean.terms = mterms, dispersion.terms = dterms, 
                       xlevels = xlevels, zlevels = zlevels, mean.contrasts = attr(x, "contrasts"), 
                       dispersion.contrasts = attr(z, "contrasts"), na.action = nmissing, 
                       init.mean.fit = init.mean.fit, init.dispersion.fit = init.dispersion.fit, 
                       matrix.inverse.method = minv.method))
    class(fit) <- c("jmdem")
    return(fit)

}

#*********************************************************************#
#**************** jmdem.prefit (JMDEM pre-estimation) ****************#
#*********************************************************************#
jmdem.prefit <- function(y, x, z, weights, start = NULL, etastart = NULL, mustart = NULL, phistart = NULL, 
                         moffset = NULL, doffset = NULL, mfamily, dfamily, control = list(), mintercept = TRUE,
                         dintercept = TRUE, nmterms, ndterms, xnames, znames, dev.type = c("deviance", "pearson"),
                         betastart = NULL, lambdastart = NULL, beta.first = TRUE, na.action, ...)
{
    ## check for multiple objects
    orgmatch <- match.call()
    arglist <- c(as.list(orgmatch)[-1], list(...))
    argnam <- names(arglist)

    ### Check whether some character inputs are correct
    dev.type <- match.arg(dev.type)
    
    ### Create the function argument lists of glm.fit
    arg.glm.fit <- as.list(args(glm.fit))
    arg.glm.fit <- arg.glm.fit[!names(arg.glm.fit) == "" & !names(arg.glm.fit) == "..."]

    ### Sample size
    n <- nrow(as.matrix(y))

    ### Get the glm.control arguments
    null.approx <- if (!is.null(arglist$control$null.approx)) {arglist$control$null.approx} else {1e-08}
    arglist$control <- arglist$control[!names(arglist$control) == "trace"]
    names(arglist$control)[names(arglist$control) == "prefit.trace"] <- "trace"
    arglist$control <- if (is.null(arglist$control)) {glm.control()} else {arglist$control <- arglist$control[which(names(arglist$control) %in% c("epsilon", "maxit", "trace"))]}

    ### Set flag parameters
    binit.done <- FALSE
    linit.done <- FALSE
    
    ### Estimate initial parameter values
    for (j in 1L:2L)
    {
        ### Estimate dispersion model parameters
        if ((beta.first & j == 1) | (!beta.first & j == 2))
        {
            ### Estimate initial value of beta
            if (is.null(betastart))
            {
                ### Get the GLM arguments
                marglist <- c(lwhich(arglist, arg.glm.fit), list(family = mfamily, offset = moffset, intercept = mintercept))

                # Estimation
                tempmodel <- try(do.call(glm.fit, args = marglist), silent = TRUE)
                if (length(which(class(tempmodel) %in% "try-error")))
                {
                    beta <- setNames(rep(0, nmterms), xnames)
                    eta <- x %*% beta + moffset
                    mu <- mfamily$linkinv(eta)
                } 
                else 
                {
                    beta <- coef(tempmodel)
                    beta[is.na(beta)] <- 0
                    eta <- tempmodel$linear.predictors
                    mu <- fitted(tempmodel)
                }
            } 
            else 
            {
                # Transfer predefined parameter values
                beta <- setNames(betastart, xnames)
                eta <- x %*% beta + moffset
                mu <- mfamily$linkinv(eta)
            }
            binit.done <- TRUE
            if (binit.done & linit.done) {break}
        }

        ### Estimate dispersion model parameters
        if ((beta.first & j == 2) | (!beta.first & j == 1))
        {
            ### Compute temporary deviance to estimate of initial value of lambda
            if (j == 1) {mu <- arglist$mustart}
            if (dev.type == "deviance") 
            {
                tempdev <- mfamily$dev.resids(y = y, mu = mu, wt = arglist$weights)
            }
            else if (dev.type == "pearson") 
            {
                tempdev <- arglist$weights * (y - mu) ^ 2 / mfamily$variance(mu)
            }
            tempdev[tempdev == 0] <- null.approx
            
            ### Estimate initial value of lambda
            if (is.null(lambdastart)) 
            {
                ### Get the GLM arguments
                darglist <- list(y = tempdev, x = z, family = arglist$dfamily, weights = rep(0.5, n), offset = doffset, intercept = dintercept)

                # Estimation
                tempmodel <- try(do.call(glm.fit, args = darglist), silent = TRUE)
                if (length(which(class(tempmodel) %in% "try-error"))) 
                {
                    lambda <- setNames(rep(1, ndterms), znames)
                    delta <- z %*% lambda + doffset
                    phi <- dfamily$linkinv(delta)
                } 
                else 
                {
                    lambda <- coef(tempmodel)
                    lambda[is.na(lambda)] <- 0
                    delta <- tempmodel$linear.predictors
                    phi <- fitted(tempmodel)
                }
            } 
            else 
            {
                
                # Transfer predefined parameter values
                lambda <- setNames(lambdastart, znames)
                delta <- z %*% lambda + doffset
                phi <- dfamily$linkinv(delta)
            }
            linit.done <- TRUE
            if (binit.done & linit.done) {break}
        }
    
    }
    return(list(beta = beta, eta = eta, mu = mu, lambda = lambda, delta = delta, phi = phi))
}
    
#*********************************************************************#
#******************** jmdem.fit (JMDEM estimation) *******************#
#*********************************************************************#
jmdem.fit <- function(x, y, z = NULL, weights, mfamily = gaussian, dfamily = Gamma, 
                      mu, phi, beta, lambda, moffset = NULL, doffset = NULL, 
                      dev.type = c("deviance", "pearson"), hessian = TRUE, 
                      method = c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN", "Brent"), 
                      disp.adj = FALSE, df.adj = FALSE, full.loglik = FALSE, 
                      control = list(), mintercept = TRUE, dintercept = TRUE, 
                      grad.func = TRUE, lower = -Inf, upper = Inf, ...)

{
    ### Listing of all model arguments
    orgmatch <- match.call()
    arglist <- as.list(orgmatch)[-1]
    argnam <- names(arglist)

    ### Check whether some character inputs are correct
    dev.type <- match.arg(dev.type)
    method <- match.arg(method)
    
    ### Control list adjustment
    est.tol <- arglist$control$epsilon
    est.maxit <- arglist$control$maxit
    est.trace <- arglist$control$fit.trace
    null.approx <- arglist$control$null.approx <- if (is.null(arglist$control$null.approx)) {1e-08} else {arglist$control$null.approx}
    arglist$control <- arglist$control[names(arglist$control)[!(names(arglist$control) %in% c("epsilon", "prefit.trace", "fit.trace", "null.approx"))]]

    ### Variable names
    ynam <- names(y)
    xvarnam <- colnames(x)
    zvarnam <- colnames(z)
    
    ### Sample size
    n <- nrow(as.matrix(y))

    ### Set empty flag if no independent variables are specified
    nvars <- ncol(x) + ncol(z)
    empty <- nvars == 0

    ### Initial parameter value
    par <- c(beta, lambda)
    npar <- length(par)
    nbeta <- length(beta)
    nlambda <- length(lambda)
    beta.index <- if (nbeta) {1L:nbeta} else {numeric(0)}
    lambda.index <- if (nlambda) {(nbeta + 1):npar} else {numeric(0)}

    ### QR decomposition of the model matrices
    mqr <- qr(x)
    dqr <- qr(z)
    mR <- qr.R(mqr)
    dR <- qr.R(dqr)
    mqr$qr <- as.matrix(mqr$qr)
    dqr$qr <- as.matrix(dqr$qr)
    
    ### Rank of the model matrices
    mrank <- mqr$rank
    drank <- dqr$rank

    ### Effects of the model matrices
    meffects <- setNames(qr.qty(mqr, y), c(names(beta), rep("", n - nbeta)))
    deffects <- setNames(qr.qty(dqr, y), c(names(lambda), rep("", n - nlambda)))

    ### Unify the specification of the deviance type
    dev.type <- tolower(dev.type)
    
    ### Save the prior weights
    prior.weights <- arglist$weights
    
    ### Check and extend call list by offset
    if (is.null(moffset)) {moffset <- rep(0L, n)}
    if (is.null(doffset)) {doffset <- rep(0L, n)}
    
    ### Create the function argument lists of loglik
    arg.loglik <- as.list(args(loglik)) #Function loglik
    arg.loglik <- arg.loglik[!names(arg.loglik) == "" & !names(arg.loglik) == "..."]

    ### Create the function argument lists of optim
    arg.optim <- as.list(args(optim)) #Function optim
    arg.optim <- arg.optim[!names(arg.optim) == "" & !names(arg.optim) == "..."]

    ### Create the user argument lists
    larglist <- lwhich(arglist, arg.loglik) #loglik function arguments 
    oarglist <- lwhich(arglist, arg.optim) #original arguments
    uarglist <- lwhich(arg.optim, arglist, reverse = TRUE) #user arguments

    ### Transfer vector indices for beta and gamma
    larglist$beta.index <- beta.index
    larglist$lambda.index <- lambda.index
    
    ### Define function name for optim
    if (length(which(names(arglist) %in% "fn")) == 0L){uarglist$fn <- quote(loglik)}

    ### Evaluate control parameters for optim
    if (!length(which(names(arglist) %in% "control")) == 0L) {oarglist$control <- lapply(oarglist$control, "eval", envir = environment())}
    
    ### If optimised by "L-BFGS-B"
    if (method == "L-BFGS-B") 
    {
        # Delete "reltol" and "abstol" from the list 
        oarglist$control <- lwhich(oarglist$control, oarglist$control[c("reltol", "abstol")], reverse = TRUE)
        
        # Delete "reltol" and "abstol" from the list 
        uarglist$lower <- rep(as.numeric(lower), npar)
        uarglist$upper <- rep(as.numeric(upper), npar)
    }
    ### Finalise argument list for estimation
    farglist <- append(larglist, append(uarglist, oarglist))

    ### Add gradient function is required
    farglist$null.approx <- null.approx
    
    ### Add gradient function is required
    if (grad.func) {farglist$gr <- gloglik}

    ### Finalise argument list for estimation
    dev.old <- n * npar * 1000

    ### Set flag for convergence
    converged <- FALSE

    ### Estimation
    for (i in 1L:est.maxit)
    {
        ### Provide the parameter values for optimisation
        farglist$par <- par

        ### Update the nested variable
        if (length(which("eta" %in% colnames(z))))
        {
            z[, "eta"] <- x %*% beta + moffset
            farglist$z <- z
        }
        if (length(which("delta" %in% colnames(x))))
        {
            x[, "delta"] <- z %*% lambda + doffset
            farglist$x <- x
        }
        
        ### Fitting model
        fit <- do.call("optim", farglist)

        ### Get model estimates
        par <- fit$par
        beta <- if (nbeta) {par[beta.index]} else {numeric(0)}
        lambda <- if (nlambda) {par[lambda.index]} else {numeric(0)}

        ### Get mean and dispersion effects
        eta <- if (nbeta) {setNames(as.numeric(x %*% beta + moffset), ynam)} else {rep(0L, n)}
        delta <- if (nlambda) {setNames(as.numeric(z %*% lambda + doffset), ynam)} else {rep(0L, n)}
        mu <- setNames(as.numeric(mfamily$linkinv(eta)), ynam)
        phi <- setNames(as.numeric(dfamily$linkinv(delta)), ynam)
        comp <- loglik.comp(family = mfamily, y = y, mu = mu, weights = prior.weights, phi = phi)
        if (disp.adj & dev.type == "pearson" & any(comp$cumadj < 0)) {stop("Dispersion adjustment (disp.adj) leads to negative working weights. Please rerun the function with disp.adj = FALSE.")}
        
        ### Dispersion adjustment
        farglist$weights <- work.weights <- prior.weights / (if (disp.adj & dev.type == "pearson") {comp$cumadj} else {1})

        ### Compute re-weights in the estimation
        weights <- sqrt(work.weights * (mfamily$mu.eta(eta) ^ 2) * (1 / mfamily$variance(mu)))
        
        ### Compute the updated deviance
        deviance.residuals <- setNames(as.numeric(sqrt(mfamily$dev.resids(y = y, mu = mu, wt = work.weights)) * sign(y - mu)), ynam)
        pearson.residuals <- setNames(as.numeric((y - mu) * sqrt(work.weights / mfamily$variance(mu))), ynam)
        residuals <- if (dev.type == "deviance") {deviance.residuals} else if (dev.type == "pearson") {residuals <- pearson.residuals}
        dev.new <- sum(residuals ^ 2)

        ### Show fit trace if required
        if (est.trace) {cat("Deviance =", round(dev.new, 4), "Iterations - ", i, "\n", sep = " ")}
        
        ### Interrupt the estimation algorithm if converge
        dev.diff <- abs(dev.new - dev.old) / (0.1 + abs(dev.new))
        if (dev.diff < est.tol)
        {
            converged <- TRUE
            break
        }
        dev.old <- dev.new
    }
    
    ### Show warning if not converged
    if (!converged) {warning("jmdem.fit: algorithm did not converge", call. = FALSE)}

    ### Degree of freedoms and rank of the fitted model
    n.ok <- n - sum(weights == 0)
    rank <- if (empty) {0} else {mrank + drank}
    resdf <- max(0, n.ok - rank)

    ### Calculate individual likelihood
    larglist <- lwhich(farglist, arg.loglik)
    iloglik <- do.call("loglik", append(larglist, list(allobs = TRUE)))
    if (length(par)) {iloglik <- setNames(as.numeric(iloglik), ynam)}

    ### AIC of the fitted model
    aic <- - 2 * fit$value + 2 * rank
    
    ### Include hessian matrix
    imat <- if (hessian) {-fit$hessian} else {do.call("fisherinfo", larglist)}

    ### Create return list
    res <- list(coefficients = par, beta = beta, lambda = lambda, residuals = residuals, deviance.residuals = deviance.residuals, 
                pearson.residuals = pearson.residuals, fitted.values = mu, dispersion = phi, mean.R = mR, dispersion.R = dR, 
                mean.effects = meffects, dispersion.effects = deffects, mean.rank = mrank, dispersion.rank = drank, rank = rank, 
                mean.qr = structure(mqr[c("qr", "rank", "qraux", "pivot")], class = "qr"), 
                dispersion.qr = structure(dqr[c("qr", "rank", "qraux", "pivot")], class = "qr"), 
                mean.family = mfamily, dispersion.family = dfamily, mean.linear.predictors = eta, dispersion.linear.predictors = delta, 
                deviance = dev.new, individual.loglik = iloglik, aic = aic, iter = i, weights = weights, 
                prior.weights = work.weights, info.matrix = imat, 
                df.residual = resdf, y = y, x = x, z = z, log.llh = fit$value, 
                converged = converged, gradient = grad.func, deviance.type = dev.type, 
                information.type = if (hessian) {"Hessian"} else {"Fisher"}, dispersion.adjustment = disp.adj, 
                df.adjustment = df.adj, optim.method = method)
    return(res)
}

#*********************************************************************#
#************ residual.jmdem (compute residuals of JMDEM) ************#
#*********************************************************************#
residuals.jmdem <- function(object, type = c("deviance", "pearson", "working",
                                             "response", "partial"), ...)
{
    ### Check whether some character inputs are correct
    type <- match.arg(type)

    ### Get information from object
    y <- object$y
    ynam <- names(y)
    r <- object$residuals
    mu <- object$fitted.values
    wts <- object$prior.weights
    switch(type,
           deviance=,pearson=,response=
               if(is.null(y)) {
                   mu.eta <- object$family$mu.eta
                   eta <- object$linear.predictors
                   y <-  mu + r * mu.eta(eta)
               })
    
    ### Compute residuals
    res <- switch(type,
                  deviance = if(object$df.residual > 0) {
                      d.res <- setNames(as.numeric(sqrt(pmax((object$mean.family$dev.resids)(y, mu, wts), 0)) * sign(y - mu)), ynam)
                  } else rep.int(0, length(mu)),
                  pearson = setNames(as.numeric((y - mu) * sqrt(wts / object$mean.family$variance(mu))), ynam),
                  working = r,
                  response = y - mu,
                  partial = r
    )

    if (!is.null(object$na.action)) {res <- naresid(object$na.action, res)}
    if (type == "partial") {res <- res + predict(object, type = "terms")} ## need to avoid doing naresid() twice.
        
    return(res)
    
}

#*********************************************************************#
#* family.jmdem (extract the distribution family of a JMDEM response) #
#*********************************************************************#
family.jmdem <- function(object, submodel = c("both", "mean", "dispersion"), ...) 
{
    ### Check whether some character inputs are correct
    submodel <- match.arg(submodel)

    ### Extract families of different submodels
    if (submodel == "mean")
    {
        res <- object$mean.family
    }
    else if (submodel == "dispersion")
    {
        res <- object$dispersion.family
    }
    else if (submodel == "both")
    {
        res <- list(mean.family = object$mean.family, dispersion.family = object$dispersion.family)
    }
    return(res)
}

#*********************************************************************#
#****** formula.jmdem (extract the formulas of a JMDEM object) *******#
#*********************************************************************#
formula.jmdem <- function(x, submodel = c("both", "mean", "dispersion"), ...)
{
    ### Check whether some character inputs are correct
    submodel <- match.arg(submodel)

    ### Extract formulas of different submodels
    if (submodel == "mean" | submodel == "both")
    {
        if (!is.null(x$mean.formula))
        {
            m.form <- x$mean.terms
            environment(m.form) <- environment(x$mean.formula)
        }
        else
        {
            m.form <- formula(x$mean.terms)
        }
    }
    if (submodel == "dispersion" | submodel == "both")
    {
        if (!is.null(x$dispersion.formula))
        {
            d.form <- x$dispersion.terms
            environment(d.form) <- environment(x$dispersion.formula)
        }
        else
        {
            d.form <- formula(x$dispersion.terms)
        }
    }
    return(if (submodel == "mean") {m.form} else if (submodel == "dispersion") {d.form} else if (submodel == "both") {list(mean.formula = m.form, dispersion.formula = d.form)})
}

#*********************************************************************#
#********* model.frame.jmdem (extract design matrix from fit) ********#
#*********************************************************************#
model.matrix.jmdem <- function (object, submodel = c("both", "mean", "dispersion"), ...)
{
    ### Check whether some character inputs are correct
    submodel <- match.arg(submodel)
    
    ### Extract design matrix
    if (submodel == "mean" || submodel == "both")
    {
        m.res <- model.matrix(object = object$mean.formula, data = object$mean.model)
    }
    if (submodel == "dispersion" || submodel == "both")
    {
        d.res <- model.matrix(object = object$dispersion.formula, data = object$dispersion.model)
    }
    
    return(if (submodel == "mean") {m.res} else if (submodel == "dispersion") {d.res} else if (submodel == "both") {list(mean = m.res, dispersion = d.res)})    
}

#*********************************************************************#
#*********** model.frame.jmdem (extract data frame for fit) **********#
#*********************************************************************#
model.frame.jmdem <- function (formula, submodel = c("both", "mean", "dispersion"), ...)
{
    ### Check whether some character inputs are correct
    submodel <- match.arg(submodel)

    ### Listing of all function arguments
    orgmatch <- match.call()
    arglist <- as.list(orgmatch)[-1]
    argnam <- names(arglist)

    ### Select only the relevant options for this function
    nargs <- arglist[match(c("data", "na.action", "subset"), argnam, 0L)]
    newdat <- as.logical(match("data", argnam, 0L))
    dnew <- if (newdat) {eval(arglist$data, parent.frame())} else {NULL}

    ### Extract the model.frame from a data set if formula is not a jmdem object
    if (class(formula) == "formula") 
    {
        if (match("delta", attr(terms.formula(formula), "term.labels"), 0L)) {dnew <- cbind(dnew, delta = rep(0, nrow(dnew)))}
        if (match("eta", attr(terms.formula(formula), "term.labels"), 0L)) {dnew <- cbind(dnew, eta = rep(0, nrow(dnew)))}
        arglist$data <- dnew
        return(do.call("model.frame", arglist))
    }
    
    ### Extract model frames of different submodels
    if (submodel == "mean" || submodel == "both")
    {
        if (length(nargs) || is.null(formula$mean.model)) 
        {
            fcall <- formula$call
            fcall$fit.method <- "model.frame"
            fcall[[1L]] <- quote(jmdem)
            fcall[names(nargs)] <- nargs
            if (match("delta", attr(formula$mean.terms, "term.labels"), 0L)) {fcall$"mformula" <- update(formula$mean.formula, . ~ . - delta)}
            env <- environment(formula$mean.terms)
            if (is.null(env)) {env <- parent.frame()}
            temp.res <- eval(fcall, env)
            m.res <- temp.res$mean
            if (match("delta", attr(formula$mean.terms, "term.labels"), 0L)) 
            {
                z <- model.matrix(formula$dispersion.formula, data = if (newdat) {dnew} else {temp.res$dispersion})
                m.res$delta <- z %*% as.matrix(formula$lambda) + if (is.null(model.offset(temp.res$dispersion))) {rep(0, nrow(z))} else {model.offset(temp.res$dispersion)}
            }
        }
        else
        {
            m.res <- formula$mean.model
        }
    }
    if (submodel == "dispersion" || submodel == "both")
    {
        if (length(nargs) || is.null(formula$dispersion.model)) 
        {
            fcall <- formula$call
            fcall$fit.method <- "model.frame"
            fcall[[1L]] <- quote(jmdem)
            fcall[names(nargs)] <- nargs
            if (match("eta", attr(formula$dispersion.terms, "term.labels"), 0L)) {fcall$"dformula" <- update(formula$dispersion.formula, ~ . - eta)}
            env <- environment(formula$dispersion.terms)
            if (is.null(env)) {env <- parent.frame()}
            temp.res <- eval(fcall, env)
            d.res <- temp.res$dispersion
            if (match("eta", attr(formula$dispersion.terms, "term.labels"), 0L)) 
            {
                x <- model.matrix(formula$mean.formula, data = if (newdat) {dnew} else {temp.res$mean})
                d.res$eta <- x %*% as.matrix(formula$beta) + if (is.null(model.offset(temp.res$mean))) {rep(0, nrow(x))} else {model.offset(temp.res$mean)}
            }
        }
        else
        {
            d.res <- formula$dispersion.model
        }
    }

    return(if (submodel == "mean") {m.res} else if (submodel == "dispersion") {d.res} else if (submodel == "both") {list(mean = m.res, dispersion = d.res)})
}

#*********************************************************************#
#************** update.frame.jmdem (update a JMDEM fit) **************#
#*********************************************************************#
update.jmdem <- function (object, mformula, dformula, ...)
{
    
    ### Listing of all function arguments
    orgmatch <- match.call()
    arglist <- as.list(orgmatch)[-1]
    arglist <- arglist[!names(arglist) == "object"]
    
    ### Update the mean and dispersion formulas
    arglist$mformula <- update(object$mean.formula, if (match("mformula", names(arglist), 0L)) {mformula} else {. ~ .})
    arglist$dformula <- update(object$dispersion.formula, if (match("dformula", names(arglist), 0L)) {dformula} else { ~ .})
    
    ### Update model
    fcall <- object$call
    fcall[[1L]] <- quote(jmdem)
    fcall[names(arglist)] <- arglist
    env <- environment(object$mean.terms)
    if (is.null(env)) {env <- parent.frame()}
    
    res <- eval(fcall, env)
    return(res)
    
}

#*********************************************************************#
#**************** predict.jmdem (prediction by JMDEM) ****************#
#*********************************************************************#
predict.jmdem <- function(object, newdata = NULL, type = c("link", "response"), 
                          se.fit = FALSE, na.action = na.pass, ...) 
{
    type <- match.arg(type)

    na.act <- object$na.action
    object$na.action <- NULL
    if (!se.fit)
    {
        if (missing(newdata)) 
        {
            pred <- switch(type, link = object$mean.linear.predictors, response = object$fitted.values)
            if (!is.null(na.act)) {pred <- napredict(na.act, pred)}
        }
        else 
        {
            pred <- newdata.predict.jmdem(object, newdata, se.fit, type = ifelse(type == "link", "response", type), na.action = na.action)
            switch(type, response = {pred <- family(object)$mean.family$linkinv(pred)}, link = )
        }
    }
    else 
    {
        pred <- newdata.predict.jmdem(object, newdata, se.fit, type = ifelse(type == "link", "response", type), na.action = na.action)
        fit <- pred$fit
        se.fit <- pred$se.fit
        switch(type, response = {se.fit <- se.fit * abs(family(object)$mean.family$mu.eta(fit))
                                 fit <- family(object)$mean.family$linkinv(fit)}, link = )
        if (missing(newdata) && !is.null(na.act)) 
        {
            fit <- napredict(na.act, fit)
            se.fit <- napredict(na.act, se.fit)
        }
        pred <- list(fit = fit, se.fit = se.fit)
    }
    
    return(pred)
}

#*********************************************************************#
#**** newdata.predict.jmdem (prediction code for new observations) ***#
#*********************************************************************#
newdata.predict.jmdem <- function (object, newdata, se.fit = FALSE, scale = NULL, df = Inf, 
                           type = c("response"), na.action = na.pass, ...) 
{
    if (!inherits(object, "jmdem")) {warning("calling predict.lm(<fake-lm-object>) ...")}
    mt <- terms(object$mean.model)
    dt <- terms(object$dispersion.model)
    n <- length(object$residuals)
    
    if (missing(newdata) || is.null(newdata)) 
    {
        mm <- X <- model.matrix(object)$mean
        md <- Z <- model.matrix(object)$dispersion
        mmDone <- TRUE
        m.offset <- object$mean.offset
        d.offset <- object$dispersion.offset
    }
    else 
    {
        mterms <- delete.response(mt)
        dterms <- delete.response(dt)
        if (match("delta", attr(mterms, "term.labels"), 0L)) {newdata$delta <- rep(0, nrow(newdata))}
        if (match("eta", attr(dterms, "term.labels"), 0L)) {newdata$eta <- rep(0, nrow(newdata))}
        m <- model.frame(mterms, newdata, na.action = na.action, xlev = object$xlevels)
        d <- model.frame(dterms, newdata, na.action = na.action, xlev = object$zlevels)

        if (!is.null(cl <- attr(mterms, "dataClasses"))) {.checkMFClasses(cl, m)}
        X <- model.matrix(mterms, m, contrasts.arg = object$mean.contrasts)
        m.offset <- rep(0, nrow(X))
        if (!is.null(off.num <- attr(mt, "offset"))) 
        {
            for (i in off.num) 
            {
                m.offset <- m.offset + eval(attr(mt, "variables")[[i + 1]], newdata)
            }
        }
        if (!is.null(object$call$moffset)) {m.offset <- m.offset + eval(object$call$moffset, newdata)}

        if (!is.null(cl <- attr(dterms, "dataClasses"))) .checkMFClasses(cl, d)
        Z <- model.matrix(dterms, d, contrasts.arg = object$dispersion.contrasts)
        d.offset <- rep(0, nrow(Z))
        if (!is.null(doff.num <- attr(dt, "offset"))) 
        {
            for (i in doff.num) 
            {
                d.offset <- d.offset + eval(attr(dt, "variables")[[i + 1]], newdata)
            }
        }
        if (!is.null(object$call$doffset)) {d.offset <- d.offset + eval(object$call$doffset, newdata)}
        mmDone <- FALSE
        
    }
    
    mp <- object$mean.rank
    mp1 <- seq_len(mp)
    mpiv <- if (mp) {object$mean.qr$pivot[mp1]}
    if (mp < ncol(X) && !(missing(newdata) || is.null(newdata))) {warning("prediction from a rank-deficient fit may be misleading")}
    beta <- object$beta
    m.predictor <- drop(X[, mpiv, drop = FALSE] %*% beta[mpiv] + if (!is.null(m.offset)) {m.offset} else {rep(0, nrow(X))})

    dp <- object$dispersion.rank
    dp1 <- seq_len(dp)
    dpiv <- if (dp) {object$dispersion.qr$pivot[dp1]}
    if (dp < ncol(Z) && !(missing(newdata) || is.null(newdata))) {warning("prediction from a rank-deficient fit may be misleading")}
    lambda <- object$lambda
    d.predictor <- drop(Z[, dpiv, drop = FALSE] %*% lambda[dpiv] + if (!is.null(d.offset)) {d.offset} else {rep(0, nrow(Z))})

    if (!(missing(newdata) || is.null(newdata)) )
    {
        if (match("delta", attr(mt, "term.labels"), 0L)) 
        {
            newdata$delta <- d.predictor
            X[, "delta"] <- d.predictor
            m.predictor <- drop(X[, mpiv, drop = FALSE] %*% beta[mpiv] + if (!is.null(m.offset)) {m.offset} else {rep(0, nrow(X))})
        }
        if (match("eta", attr(dt, "term.labels"), 0L)) 
        {
            newdata$eta <- m.predictor
            Z[, "eta"] <- m.predictor
            d.predictor <- drop(Z[, dpiv, drop = FALSE] %*% lambda[dpiv] + if (!is.null(d.offset)) {d.offset} else {rep(0, nrow(Z))})
        }
    }
    
    type <- match.arg(type)
    cmiv <- c(mpiv, dpiv + max(object$mean.qr$pivot))
    if (se.fit) 
    {
        im <- object$info.matrix
        covmat <- minv(im, object$matrix.inverse.method)$inv
        XRinv <- cbind(X[, mpiv, drop = FALSE], Z[, dpiv, drop = FALSE])
        predict_var <- diag(XRinv %*% covmat[cmiv, cmiv] %*% t(XRinv))
        se <- sqrt(predict_var)
    }

    if (missing(newdata) && !is.null(na.act <- object$na.action)) 
    {
        m.predictor <- napredict(na.act, m.predictor)
        if (se.fit) {se <- napredict(na.act, se)}
    }

    return(if (se.fit) {list(fit = m.predictor, se.fit = se)} else {m.predictor})
}

#*********************************************************************#
#************ jmdem.control (control parameters for JMDEM) ***********#
#*********************************************************************#
jmdem.control <- function(maxit = 100, epsilon = 1e-8, prefit.trace = FALSE, 
                          fit.trace = FALSE, null.approx = 1e-8, trace = 0, 
                          fnscale = -1, parscale = 1, ndeps = 0.001, 
                          abstol = -Inf, reltol = sqrt(.Machine$double.eps), 
                          alpha = 1, beta = 0.5, gamma = 2, REPORT = 10, 
                          type = 1, lmm = 5, factr = 1e+07, pgtol = 0, 
                          temp = 10, tmax = 10)
{
    ### Listing of all control arguments
    orgmatch <- match.call()
    arglist <- as.list(orgmatch)[-1]
    
    ### Listing of all model arguments
    conlist <- as.list(args(jmdem.control))[!names(as.list(args(jmdem.control))) == ""]
    
    ### Mix the user and default arguments
    if (length(arglist) > 0L)
    {
        ans <- append(arglist, lwhich(conlist, arglist, reverse = TRUE))
    } 
    else 
    {
        ans <- conlist
    }
    return(ans)
}

#*********************************************************************#
#*************** print.jmdem (print short JMDEM output) **************#
#*********************************************************************#
print.jmdem <- function(x, digits = max(3L, getOption("digits") - 3L), ...)
{
    ### Get scaled deviance
    scaled.null.deviance <- sum(x$residuals.null ^ 2 / x$dispersion.null)
    scaled.deviance <- sum(x$residuals ^ 2 / x$dispersion)
    
    ### Get significant digits
    dig.beta <- digits + max(floor(log10(abs(x$beta))) + 1)
    dig.lambda <- digits + max(floor(log10(abs(x$lambda))) + 1)
    dig.dev <- digits + max(floor(log10(scaled.deviance)) + 1, floor(log10(scaled.null.deviance)) + 1)
    dig.aic <- digits + floor(log10(abs(x$aic))) + 1

    ### Print "Call:"
    cat("\nCall:  ", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
    
    ### Print "Coefficients:"
    if (length(coef(x))) 
    {
        cat("Mean Effect Coefficients")
        if (is.character(co <- x$mean.contrasts))
        {
            cat("  [contrasts: ", apply(cbind(names(co), co), 1L, paste, collapse = "="), "]")
        }
        cat(":\n")
        print.default(format(x$beta, digits = dig.beta), print.gap = 2, quote = FALSE)
        cat("\nDispersion Effect Coefficients")
        if (is.character(co <- x$mean.contrasts))
        {
            cat("  [contrasts: ", apply(cbind(names(co), co), 1L, paste, collapse = "="), "]")
        }
        cat(":\n")
        print.default(format(x$lambda, digits = dig.lambda), print.gap = 2, quote = FALSE)
    } 
    else 
    {
        cat("No coefficients\n\n")
    }
    
    ### Print "Degrees of Freedom:"
    cat("\nDegrees of Freedom:", x$df.null, "Total (i.e. Null); ", x$df.residual, "Residual\n")

    ### Print na.action
    if (nzchar(mess <- naprint(x$na.action))) {cat("  (",mess, ")\n", sep = "")}
    
    ### Print deviance and AIC
    cat("Scaled Null Deviance:	    ", format(signif(scaled.null.deviance, dig.dev)), "\nScaled Residual Deviance: ", 
        format(signif(scaled.deviance, dig.dev)), "\tAIC:", format(signif(x$aic, dig.aic)))
    cat("\n")
    
    invisible(x)
    
}

#*********************************************************************#
#******************* summary.jmdem (JMDEM summary) *******************#
#*********************************************************************#
summary.jmdem <- function(object, correlation = FALSE, symbolic.cor = FALSE, ...)
{
    ### Listing of all function arguments
    orgmatch <- match.call()
    arglist <- as.list(orgmatch)[-1]

    ### Set default digits and scientific notation settings
    if (!is.null(arglist$digits)) {digits = arglist$digits} else {digits <- max(3L, getOption("digits") - 3L)}
    if (!is.null(arglist$scientific)) {scientific = arglist$scientific} else {scientific <- FALSE}
    
    ### Degrees of freedom of the residual deviance
    df.r <- object$df.residual
    
    ### Checked if any coefficients were not estimated
    aliased <- is.na(coef(object))
    
    ### Rank of the fitted model
    p <- object$rank
    mp <- object$mean.rank
    dp <- object$dispersion.rank

    ### Matrix inverse method
    minv.method <- object$matrix.inverse.method
    
    ### If the fitted model is not empty
    if (p > 0) 
    {
        ### Create the length of the valid parameters based on the rank
        mp1 <- 1L:mp
        dp1 <- 1L:dp
        p1 <- 1L:(mp + dp)
        
        ### Create the valid indices
        mpivot <- if (mp) {object$mean.qr$pivot[mp1]} else {numeric(0)}
        dpivot <- if (dp) {object$dispersion.qr$pivot[dp1]} else {numeric(0)}
        model.pivot <- c(mpivot, dpivot + length(object$beta))
        
        ### Get the valid coefficients and its covariance matrix
        coef.p <- object$coefficients[model.pivot]
        coeferr <- std.err(object, range = model.pivot, minv.method = minv.method)
        s.err <- coeferr$s.err
        covmat <- coeferr$covmat
        if (correlation) {cormat <- covmat / outer(s.err, s.err)}
        covmat.method <- coeferr$covmat.method
        tvalue <- coef.p / s.err
        dn <- c("Estimate", "Std. Error")
        pvalue <- 2 * pt(-abs(tvalue), df.r)
        
        ### Create the coefficient table
        coef.table <- cbind(coef.p, s.err, tvalue, pvalue)
        dimnames(coef.table) <- list(names(coef.p), c(dn, "t value", "Pr(>|t|)"))
        
        ### Separate the mean model coefficient tables
        if (mp)
        {
            mean.coef.table <- coef.table[mp1, , drop = FALSE]
        }
        else
        {
            mean.coef.table <- matrix(NA, 0L, 4L)
            dimnames(mean.coef.table) <- list(NULL, c("Estimate", "Std. Error", "t value", "Pr(>|t|)"))
        }
        
        ### Separate the dispersion model coefficient tables
        if (dp)
        {
            dispersion.coef.table <- coef.table[dp1 + mp, , drop = FALSE]
        }
        else
        {
            dispersion.coef.table <- matrix(NA, 0L, 4L)
            dimnames(dispersion.coef.table) <- list(NULL, c("Estimate", "Std. Error", "t value", "Pr(>|t|)"))
        }
        
        ### Full degrees of freedom
        df.f <- NCOL(object$mean.qr$qr) + NCOL(object$dispersion.qr$qr)
    } 
    else 
    {
        ### Create an empty coefficient table if no valid coefficients
        coef.table <- matrix(NA, 0L, 4L)
        mean.coef.table <- matrix(NA, 0L, 4L)
        dispersion.coef.table <- matrix(NA, 0L, 4L)
        dimnames(coef.table) <- list(NULL, c("Estimate", "Std. Error", "t value", "Pr(>|t|)"))
        dimnames(mean.coef.table) <- list(NULL, c("Estimate", "Std. Error", "t value", "Pr(>|t|)"))
        dimnames(dispersion.coef.table) <- list(NULL, c("Estimate", "Std. Error", "t value", "Pr(>|t|)"))
        covmat <- matrix(NA, 0L, 0L)
        if (correlation) {cormat <- matrix(NA, 0L, 0L)}
        df.f <- length(aliased)
        covmat.method <- NULL
    }

    ### Select the items in the jmdem object for output
    keep <- match(c("call", "mean", "dispersion", "mean.family", 
                    "dispersion.family", "deviance", "mean.terms", 
                    "dispersion.terms", "aic", "mean.contrasts", 
                    "dispersion.contrasts", "df.residual", "null.deviance", 
                    "df.null", "information.type", "iter", "na.action",
                    "dispersion.deviance", "dispersion.df.residual",
                    "dispersion.null.deviance", "dispersion.df.null",
                    "residuals.null", "dispersion.null"), 
                  names(object), 0L)
    
    ### Collect all items for output    
    ans <- c(object[keep], list(deviance.resid = object$deviance.residuals, 
                                pearson.resid = object$pearson.residuals, 
                                resid = object$residuals, coefficients = coef.table, 
                                mean.coefficients = mean.coef.table, 
                                dispersion.coefficients = dispersion.coef.table, 
                                deviance.type = object$deviance.type, aliased = aliased, 
                                df = c(object$rank, df.r, df.f), covariance = covmat, 
                                digits = digits, scientific = scientific, 
                                covmat.method = covmat.method))
    
    ### Add correlation matrix if required
    if (correlation) {
        ans$correlation <- cormat
        ans$symbolic.cor <- symbolic.cor
    }
    
    class(ans) <- "summary.jmdem"
    return(ans)
}

#*********************************************************************#
#************* print.summary.jmdem (print JMDEM summary) *************#
#*********************************************************************#
print.summary.jmdem <- function (x, digits = max(3L, getOption("digits") - 3L), 
                                 scientific = FALSE, symbolic.cor = x$symbolic.cor, 
                                 signif.stars = getOption("show.signif.stars"), ...)
{
    ### Listing of all function arguments
    orgmatch <- match.call()
    arglist <- as.list(orgmatch)[-1]
    
    ### Set default digits and scientific notation settings
    digits <- if (is.null(x$digits)) {digits} else {x$digits}
    scientific <- if (is.null(x$scientific)) {scientific} else {x$scientific}

    ### Get scaled deviance
    scaled.null.deviance <- sum(x$residuals.null ^ 2 / x$dispersion.null)
    scaled.deviance <- sum(x$resid ^ 2 / x$dispersion)
    
    ### Adjust the length of the coefficient labels
    mmaxchar <- if (nrow(x$mean.coefficients)) {max(nchar(rownames(x$mean.coefficients)))} else {0}
    dmaxchar <- if (nrow(x$dispersion.coefficients)) {max(nchar(rownames(x$dispersion.coefficients)))} else {0}
    charadd <- max(mmaxchar, dmaxchar) - c(mmaxchar, dmaxchar)
    if (nrow(x$mean.coefficients)) {rownames(x$mean.coefficients) <- paste0(rownames(x$mean.coefficients), paste(rep(" ", charadd[1]), collapse = ""), sep = "")}
    if (nrow(x$dispersion.coefficients)) {rownames(x$dispersion.coefficients) <- paste0(rownames(x$dispersion.coefficients), paste(rep(" ", charadd[2]), collapse = ""), sep = "")}

    ### Print "Call:"    
    cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
    
    ### Print "Residuals:"
    dev.type <- paste(toupper(substr(x$deviance.type, 1, 1)), tolower(substr(x$deviance.type, 2, nchar(x$deviance.type))), " Residuals: \n", sep = "")
    cat(dev.type)
    if (x$df.residual > 5) 
    {
        x$deviance.resid <- setNames(quantile(x$deviance.resid, na.rm = TRUE), c("Min", "1Q", "Median", "3Q", "Max"))
    }
    xx <- zapsmall(x$deviance.resid, digits + 1L)
    print.default(xx, digits = digits, na.print = "", print.gap = 2L)
    
    ### Print "Coefficients:"
    if (length(x$aliased) == 0L) 
    {
        cat("\nNo Coefficients\n")
    } 
    else 
    {
        df <- if ("df" %in% names(x)) x[["df"]] else NULL
        if (!is.null(df) && (nsingular <- df[3L] - df[1L]))
        {
            cat("\nCoefficients: (", nsingular, " not defined because of singularities)\n", sep = "")
        }
        
        # Mean coefficients
        cat("\nMean Effect Coefficients:\n")
        if (nrow(x$mean.coefficients))
        {
            coefs <- x$mean.coefficients
            if (!is.null(aliased <- x$aliased) && any(aliased)) 
            {
                cn <- names(aliased)
                coefs <- matrix(NA, length(aliased), 4L, dimnames = list(cn, colnames(coefs)))
                coefs[!aliased, ] <- x$mean.coefficients
            }
            printCoefmat.jmdem(coefs, digits = digits, signif.stars = TRUE, signif.legend = FALSE, na.print = "NA", scientific = scientific, ...)
        }
        else
        {
            cat("No Coefficients\n")
        }
    }
        
    ### Get significant digits of the deviances and AIC
    dig.dev1 <- if (!is.na(scaled.deviance)) {floor(log10(abs(scaled.deviance))) + 1} else {0}
    dig.dev2 <- if (!is.na(scaled.null.deviance)) {floor(log10(abs(scaled.null.deviance))) + 1} else {0}
    dig.dev <- digits + max(dig.dev1, dig.dev2)

    ### Print deviance and degrees of freedom
    cat(paste("\n", apply(cbind(paste(format(c("Scaled null", "Scaled residual"), justify = "right"), "deviance:"), 
                                format(c(scaled.null.deviance, scaled.deviance), digits = dig.dev), "on",
                                format(unlist(x[c("df.null", "df.residual")])), "degrees of freedom"), 
                          1L, paste, collapse = " "), sep = ""))
        
    if (length(x$aliased)) 
    {
        # Dispersion coefficients
        cat("\n\nDispersion Effect Coefficients:\n")
        if (nrow(x$dispersion.coefficients))
        {
            coefs <- x$dispersion.coefficients
            if (!is.null(aliased <- x$aliased) && any(aliased)) 
            {
                cn <- names(aliased)
                coefs <- matrix(NA, length(aliased), 4L, dimnames=list(cn, colnames(coefs)))
                coefs[!aliased, ] <- x$coefficients
            }
            printCoefmat.jmdem(coefs, digits = digits, signif.stars = signif.stars, na.print = "NA", scientific = scientific, ...)
        }
        else
        {
            cat("No Coefficients\n")
        }
    }

    ### Get significant digits of the deviances and AIC
    dig.dev1 <- if (!is.na(x$dispersion.deviance)) {floor(log10(abs(x$dispersion.deviance))) + 1} else {0}
    dig.dev2 <- if (!is.na(x$dispersion.null.deviance)) {floor(log10(abs(x$dispersion.null.deviance))) + 1} else {0}
    dig.dev <- digits + max(dig.dev1, dig.dev2)
    dig.aic <- digits + floor(log10(abs(x$aic))) + 1

    ### Print deviance and degrees of freedom
    cat(paste("\n", apply(cbind(paste(format(c("Scaled null", "Scaled residual"), justify = "right"), "deviance:"), 
                                format(unlist(x[c("dispersion.null.deviance", "dispersion.deviance")]), digits = dig.dev), "on",
                                format(unlist(x[c("dispersion.df.null", "dispersion.df.residual")])), "degrees of freedom"), 
                          1L, paste, collapse = " "), sep = ""))

    ### Print na.action if any missings in the data
    if (nzchar(mess <- naprint(x$na.action))) {cat("\n  (", mess, ")\n", sep = "")}

    ### Print rest information
    if (x$aic) {cat("\nAIC: ", format(x$aic, digits = dig.aic))}
    cat("\n\nNumber of iterations: ", x$iter, "\n", sep = "")
    if (!is.null(x$covmat.method)) {cat("(The covariance matrix here is the ", x$covmat.method, " of the ", if(x$information.type == "Fisher") {"Fisher Information"} else {x$information.type}, " matrix.)\n", sep = "")}

    ### Print correlation if required
    correl <- x$correlation
    if (!is.null(correl)) 
    {
        p <- NCOL(correl)
        if (p > 1) 
        {
            cat("\nCorrelation of Coefficients:\n")
            if (is.logical(symbolic.cor) && symbolic.cor) 
            {
                print(symnum(correl, abbr.colnames = NULL)) # NULL < 1.7.0 objects
            } 
            else 
            {
                correl <- format(round(correl, digits), nsmall = 2L, digits = digits)
                correl[!lower.tri(correl)] <- ""
                print(correl[-1, -p, drop=FALSE], quote = FALSE)
            }
        }
    }
    
    cat("\n")
    invisible(x)
}

#*********************************************************************#
#************** printCoefmat (print coefficient table) ***************#
#*********************************************************************#
printCoefmat.jmdem <- function(x, digits = max(3L, getOption("digits") - 2L), signif.stars = getOption("show.signif.stars"), 
                               signif.legend = signif.stars, cs.ind = 1:k, tst.ind = k + 1, zap.ind = integer(), scientific = FALSE,
                               dig.tst = max(1L, min(5L, digits) - 1), eps.Pvalue = .Machine$double.eps, na.print = "NA", 
                               P.values = NULL, has.Pvalue = nc >= 4L && substr(colnames(x)[nc], 1L, 3L) %in% c("Pr(", "p-v"), ...)
{
    ### Check whether it should proceed
    if (is.null(d <- dim(x)) || length(d) != 2L) {stop("'x' must be coefficient matrix/data frame")}
    
    ### Number of column
    nc <- d[2L]
    
    ### Stop the summary if p value options are wrong
    if (is.null(P.values)) 
    {
        scp <- getOption("show.coef.Pvalues")
        if (!is.logical(scp) || is.na(scp)) 
        {
            warning("option \"show.coef.Pvalues\" is invalid: assuming TRUE")
            scp <- TRUE
        }
        P.values <- has.Pvalue && scp
    } 
    else if (P.values && !has.Pvalue) 
    {
        stop("'P.values' is TRUE, but 'has.Pvalue' is not")
    }
    
    ### Prepare output matrix
    if (has.Pvalue && !P.values) 
    {
        d <- dim(xm <- data.matrix(x[, -nc, drop = FALSE]))
        nc <- nc - 1
        has.Pvalue <- FALSE
    } 
    else 
    {
        xm <- data.matrix(x)
        k <- nc - has.Pvalue - (if (missing(tst.ind)) {1} else {length(tst.ind)})
    }
    
    ### Stop if the wrong indices are given
    if (!missing(cs.ind) && length(cs.ind) > k) {stop("wrong k / cs.ind")}
    
    ### Initiate the final output matrix
    Cf <- array("", dim = d, dimnames = dimnames(xm))
    
    ### Check which values are non-missing
    ok <- !(ina <- is.na(xm))
    
    ### Set the cells of coefficients and std. err. rounded to zero if it is small
    for (i in zap.ind) {xm[, i] <- zapsmall(xm[, i], digits)}
    
    ### Round the decimal places of the coefficients and std. err. and truncate them to the given significant digits
    if (length(cs.ind)) 
    {
        acs <- abs(coef.se <- xm[, cs.ind, drop = FALSE])
        if (any(ia <- is.finite(acs))) 
        {
            digmin <- if (scientific) {1} else {1 + if (length(acs <- acs[ia & acs != 0])) {floor(log10(range(acs[acs != 0], finite = TRUE)))} else {0}}
            digbefdec <- max(1, digmin)
            Cf[, cs.ind] <- format(round(coef.se, max(1L, digits)), digits = digits + digbefdec, scientific = scientific)
        }
    }
    
    ### Round the decimal places of the test statistics and truncate them to the given significant digits
    if (length(tst.ind)) 
    {
        acs <- abs(tst <- xm[, tst.ind, drop = FALSE])
        if (any(ia <- is.finite(acs))) 
        {
            if (!scientific) {dig.tst <- dig.tst + 1}
            digmin <- if (scientific) {1} else {1 + if (length(acs <- acs[ia & acs != 0])) {floor(log10(range(acs[acs != 0], finite = TRUE)))} else {0}}
            digbefdec <- max(1, digmin)
            Cf[, tst.ind] <- format(round(xm[, tst.ind], digits = dig.tst), digits = dig.tst + digbefdec, scientific = scientific)
        }
    }
    
    ### Round the table if the table contains other columns
    if (any(r.ind <- !((1L:nc) %in% c(cs.ind, tst.ind, if (has.Pvalue) nc)))) {for (i in which(r.ind)) Cf[, i] <- round(xm[, i], digits = digits)}
    
    ### Get all the valid values and see whether it has to be changed to scientific if required
    ok[, tst.ind] <- FALSE
    okP <- if (has.Pvalue) {ok[, -nc]} else {ok}
    x1 <- Cf[okP]
    
    ### Use the right decimal symbol
    dec <- getOption("OutDec")
    if (dec != ".") {x1 <- chartr(dec, ".", x1)}
    
    ### Change those small values to scientific notation if required
    if (scientific) 
    {
        x0 <- (xm[okP] == 0) != (as.numeric(x1) == 0)
        if (length(not.both.0 <- which(x0 & !is.na(x0)))) {Cf[okP][not.both.0] <- format(xm[okP][not.both.0], digits = max(1L, digits - 1L))}
    }
    
    ### Format missings in the table
    if (any(ina)) {Cf[ina] <- na.print}
    
    ### Significance stars
    if (P.values) 
    {
        if (!is.logical(signif.stars) || is.na(signif.stars)) 
        {
            warning("option \"show.signif.stars\" is invalid: assuming TRUE")
            signif.stars <- TRUE
        }
    }
    
    ### Format p values with significant stars
    if (any(okP <- ok[, nc])) 
    {
        pv <- as.vector(xm[, nc])
        Cf[okP, nc] <- format.pval(round(pv[okP], 4), digits = dig.tst, eps = eps.Pvalue)
        signif.stars <- signif.stars && any(pv[okP] < 0.1)
        if (signif.stars) 
        {
            Signif <- symnum(pv, corr = FALSE, na = FALSE, cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), symbols = c("***", "**", "*", ".", " "))
            Cf <- cbind(Cf, format(Signif))
        } 
        else 
        {
            signif.stars <- FALSE
        } 
    } 
    else 
    {
        signif.stars <- FALSE
    }
    
    ### Print the table    
    print.default(Cf, quote = FALSE, right = TRUE, na.print = na.print, ...)
    
    ### Print the significant stars legend
    if (signif.stars && signif.legend) 
    {
        if ((w <- getOption("width")) < nchar(sleg <- attr(Signif, "legend"))) 
        {
            sleg <- strwrap(sleg, width = w - 2, prefix = "  ")
        }
        cat("\n---\nSignif. codes:  ", sleg, sep = "", fill = w + 4 + max(nchar(sleg, "bytes") - nchar(sleg)))
    }
    invisible(x)
}

#*********************************************************************#
#****************** anova.jmdem (ANOVA for JMDEM) ********************#
#*********************************************************************#
anova.jmdem <- function(object, ..., test = NULL, type = c("1", "3"), 
                        print.results = TRUE)
{
    ### Check whether some character input is correct
    type <- as.numeric(match.arg(arg = as.character(type), choices = c("1", "3")))
    if (!is.null(test)) {test <- match.arg(arg = test, choices = c("Rao", "Wald"))}

    ### Check for multiple objects
    dotargs <- list(...)
    named <- if (is.null(names(dotargs))) {rep_len(FALSE, length(dotargs))} else {(names(dotargs) != "")}
    if (any(named)) 
    {
        warning("the following arguments to 'anova.jmdem' are invalid and dropped: ", paste(deparse(dotargs[named]), collapse=", "))
    }
    dotargs <- dotargs[!named]
    is.jmdem <- vapply(dotargs, function(x) inherits(x, "jmdem"), NA)
    dotargs <- dotargs[is.jmdem]

    ### If more than one objects, that means compare two or more objects, then use anova.jmdemlist    
    if (length(dotargs)) {return(anova.jmdemlist(c(list(object), dotargs), test = test))}
    
    ### Score or Wald tests require a bit of extra computing
    doscore <- !is.null(test) && test == "Rao"
    dowald <- !is.null(test) && test == "Wald"
    
    ### Sample size
    n <- nrow(as.matrix(object$y))
    
    ### Prepare df and deviance columns
    resdf.x  <- numeric(0)
    resdev.x <- numeric(0)
    if (doscore || dowald) {asytest.x <- c(NA)}
    resdf.z  <- numeric(0)
    resdev.z <- numeric(0)
    if (doscore || dowald) {asytest.z <- c(NA)}
    
    ### Extract variables from model
    mvarlist <- attr(object$mean.terms, "term.labels")
    dvarlist <- attr(object$dispersion.terms, "term.labels")
    
    ### Extract formulas from model
    mformula <- as.character(object$mean.formula)
    dformula <- as.character(object$dispersion.formula)

    ### Extract response variable from model
    mresponse <- all.vars(object$mean.formula[[2]])
    
    ### Extract intercepts from model
    mintercept <- attr(object$mean.terms, "intercept")
    dintercept <- attr(object$dispersion.terms, "intercept")

    ### Determine start mean and dispersion formula
    mf <- if (type == 1) {paste(mresponse, " ~ ", mintercept, sep = "")} else if (type == 3) {paste(mformula[2], mformula[1], mformula[3], sep = " ")}
    df <- if (type == 1) {paste(" ~ ", dintercept, sep = "")} else if (type == 3) {paste(dformula[1], dformula[2], sep = " ")}

    ### Fit initial model
    mim <- mcm <- if (type == 1) {update(object, mformula = mf)} else if (type == 3) {object}
    dim <- dcm <- if (type == 1) {update(object, dformula = df)} else if (type == 3) {object}

    ### Get scaled deviance
    scaled.null.deviance <- sum(object$residuals.null ^ 2 / object$dispersion.null)
    scaled.deviance <- sum(object$residuals ^ 2 / object$dispersion)
    
    ### Initiate output table
    resdev.x <- c(resdev.x, if (type == 1) {scaled.null.deviance} else if (type == 3) {scaled.deviance})
    resdf.x <- c(resdf.x, if (type == 1) {object$df.null} else if (type == 3) {object$df.residual})
    resdev.z <- c(resdev.z, if (type == 1) {object$dispersion.null.deviance} else if (type == 3) {object$dispersion.deviance})
    resdf.z <- c(resdf.z, if (type == 1) {object$dispersion.df.null} else if (type == 3) {object$dispersion.df.residual})

    ### Create ANOVA table for mean submodel
    cf <- ""
    for (i in mvarlist)
    {
        symb <- if (type == 1) {"+"} else if (type == 3) {"-"}
        cf <- if (type == 1) {paste(cf, symb, i, sep = " ")} else if (type == 3) {paste(symb, paste(mvarlist[grep(i, mvarlist)], collapse = " - "), sep = " ")}
        uf <- as.formula(paste(mf, cf, sep = " "))
        mpm <- mcm
        mcm <- update(object, mformula = uf)
        resdev.x <- c(resdev.x, sum(mcm$residuals ^ 2 / mcm$dispersion))
        resdf.x <- c(resdf.x, mcm$df.residual)
        
        ### Run asymptotic tests if user requires
        test.list <- if (type == 1) {list(mpm, mcm)} else if (type == 3) {list(mim, mcm)} 
        asytemp <- if (doscore) {score.jmdem(test.list)} else if (dowald) {wald.jmdem(test.list)} else {NA}
        if (doscore || dowald) {asytest.x <- c(asytest.x, asytemp)}
    }
    
    ### Create ANOVA table for dispersion submodel
    cf <- ""
    for (i in dvarlist)
    {
        symb <- if (type == 1) {"+"} else if (type == 3) {"-"}
        cf <- if (type == 1) {paste(cf, symb, i, sep = " ")} else if (type == 3) {paste(symb, paste(dvarlist[grep(i, dvarlist)], collapse = " - "), sep = " ")}
        uf <- as.formula(paste(df, cf, sep = " "))
        dpm <- dcm
        dcm <- update(object, dformula = uf)
        resdev.z <- c(resdev.z, dcm$dispersion.deviance)
        resdf.z <- c(resdf.z, dcm$dispersion.df.residual)
        
        ### Run asymptotic tests if user requires
        test.list <- if (type == 1) {list(dpm, dcm)} else if (type == 3) {list(dim, dcm)} 
        asytemp <- if (doscore) {score.jmdem(test.list)} else if (dowald) {wald.jmdem(test.list)} else {NA}
        if (doscore || dowald) {asytest.z <- c(asytest.z, asytemp)}
    }
    
    ### Compute deviance
    devdf.x <- resdev.x / resdf.x
    devdf.z <- resdev.z / resdf.z
    dfdiff.x <- if (type == 1) {c(NA, -diff(resdf.x))} else if (type == 3) {c(NA, resdf.x[-1] - mim$df.residual)}
    dfdiff.z <- if (type == 1) {c(NA, -diff(resdf.z))} else if (type == 3) {c(NA, resdf.z[-1] - dim$dispersion.df.residual)}
    
    ### Contruct output tables
    table.x <- data.frame(resdf.x, resdev.x, devdf.x)
    table.z <- data.frame(resdf.z, resdev.z, devdf.z)
    if (!is.null(test)) {table.x <- cbind(table.x, dfdiff.x)}
    if (!is.null(test)) {table.z <- cbind(table.z, dfdiff.z)}

    ### Add row and column names to the tables
    tl.x <- attr(object$mean.terms, "term.labels")
    tl.z <- attr(object$dispersion.terms, "term.labels")
    if (length(tl.x) == 0L) {table.x <- table.x[1, , drop = FALSE]}
    if (length(tl.z) == 0L) {table.z <- table.z[1, , drop = FALSE]}
    dimnames(table.x) <- list(c(if (type == 1) {"NULL"} else if (type == 3) {"FULL"}, tl.x), c("Resid. Df", "Resid. Dev", "Resid. Dev / Df", if (!is.null(test)) {"Df"}))
    dimnames(table.z) <- list(c(if (type == 1) {"NULL"} else if (type == 3) {"FULL"}, tl.z), c("Resid. Df", "Resid. Dev", "Resid. Dev / Df", if (!is.null(test)) {"Df"}))
    
    ### Add asymptotic test if required
    if (doscore || dowald)
    {
        table.x <- cbind(table.x, test = asytest.x)
        table.z <- cbind(table.z, test = asytest.z)
    }
    colnames(table.x)[colnames(table.x) == "test"] <- colnames(table.z)[colnames(table.z) == "test"] <- if (doscore) {"Rao"} else if (dowald) {"Wald"}
    
    ### Create output title
    title.x <- paste0("Type ", type, " Analysis of Deviance Table\n\nResponse: ", mformula[2L], 
                      "\nMean Effect Model: ", object$mean.family$family, "(link = ", object$mean.family$link, ")",  
                      if (type == 1) {"\nTerms added sequentially (first to last)\n"} else if (type == 3) {"\nTerms removed partially\n"})
    title.z <- paste0("\nDispersion Effect Model: ", object$dispersion.family$family, "(link = ", object$dispersion.family$link, ")", 
                      if (type == 1) {"\nTerms added sequentially (first to last)\n"} else if (type == 3) {"\nTerms removed partially\n"})

    ## Calculate test statistic for F test if needed
    if (!is.null(test)) 
    {
        table.x <- stat.anova.jmdem(table = table.x, test = test)
        table.z <- stat.anova.jmdem(table = table.z, test = test)
    }
    
    ### Return the final table
    if (print.results)
    {
        print(structure(table.x, heading = title.x, digits = 4, class = c("anova", "data.frame")))
        print(structure(table.z, heading = title.z, digits = 4, class = c("anova", "data.frame")))
    }
    else
    {
        return(list(table.x = table.x, table.z = table.z))
    }    
    
}

#*********************************************************************#
#***************** anova.jmdemlist (ANOVA for JMDEM) *****************#
#*********************************************************************#
anova.jmdemlist <- function(object, ..., test = NULL)
{
    ### Determine the test to be conducted
    doscore <- !is.null(test) && test=="Rao"
    dowald <- !is.null(test) && test=="Wald"
    
    ### Find responses for all models and remove any models with a different response
    responses <- as.character(lapply(object, function(x) {deparse(x$mean.formula[[2L]])} ))
    sameresp <- responses == responses[1L]
    if (!all(sameresp)) 
    {
        object <- object[sameresp]
        warning(gettextf("models with response %s removed because response differs from model 1",
                         sQuote(deparse(responses[!sameresp]))), domain = NA)
    }
    
    ### Check whether the models are from the same dataset
    ns <- sapply(object, function(x) length(x$residuals))
    if (any(ns != ns[1L])) {stop("models were not all fitted to the same size of dataset")}
    
    ### Calculate the number of models
    nmodels <- length(object)
    if (nmodels == 1) {return(anova.jmdem(object[[1L]], test=test))}
    
    ### Extract statistics
    resdf  <- as.numeric(lapply(object, function(x) x$df.residual))
    resdev <- as.numeric(lapply(object, function(x) sum(as.numeric(x$residuals) ^ 2 / as.numeric(x$dispersion))))
    devdf <- resdev / resdf

    ### If either score or Wald test is required
    if (doscore | dowald)
    {
        # Create an empty vector with entries as the number of models
        asytest <- numeric(nmodels)
        asytest[1] <- NA
        
        # The difference of the residual degrees of freedom is calculated
        df <- -diff(resdf)
        
        # Run through all the compared models
        for (i in seq_len(nmodels-1)) 
        {
            # Test statistic with the corresponding sign
            asytest[i + 1] <- if (doscore) {score.jmdem(list(object[[i]], object[[i + 1]]))} else if (dowald) {wald.jmdem(list(object[[i]], object[[i + 1]]))}
            if (df < 0) asytest[i+1] <- -asytest[i+1]
        }
    }
    
    ### Construct table and title
    table <- data.frame(resdf, resdev, devdf, c(NA, -diff(resdf)), c(NA, -diff(resdev)))
    mvariables <- lapply(object, function(x) paste(deparse(x$mean.formula), collapse = "\n"))
    dvariables <- lapply(object, function(x) paste(deparse(x$dispersion.formula), collapse = "\n"))
    dimnames(table) <- list(1L:nmodels, c("Resid. Df", "Resid. Dev", "Resid. Dev / Df", "Df", "Deviance"))
    if (doscore) {table <- cbind(table, Rao = asytest)}
    if (dowald) {table <- cbind(table, Wald = asytest)}
    title <- "Analysis of Deviance Table\n"
    topnote <- paste(paste("Model ", format(1L:nmodels),": (Mean) ", mvariables, "; (Dispersion) ", dvariables, sep = "", collapse = "\n"), "\n", sep = "")

    ### Calculate test statistic for F test if needed
    if (!is.null(test)) 
    {
        table <- stat.anova.jmdem(table = table, test = test)
    }
    structure(table, heading = c(title, topnote), class = c("anova", "data.frame"))
}

#*********************************************************************#
#**** stat.anova.jmdem (create ANOVA table with Wald test added) *****#
#*********************************************************************#
stat.anova.jmdem <- function (table, test = c("Rao", "Wald")) 
{
    ### Check whether some character inputs are correct
    test <- if (length(test) > 1) {colnames(table)[which(colnames(table) %in% c("Rao", "Wald"))]} else {match.arg(test)}
    
    ### Assign column name according to test name
    if (length(test) == 0) {stop("Test has to be either 'Rao' or 'Wald'!")}
    dev.col <- match(test, colnames(table))
    
    ### Extend table with columns of df, test statistics, p-values
    if (sum(match(c("Rao", "Wald"), colnames(table), 0L)))
    {
        switch(test, Rao = , Wald = 
        {
            dfs <- table[, "Df"]
            vals <- table[, dev.col] * sign(dfs)
            vals[dfs %in% 0] <- NA
            vals[!is.na(vals) & vals < 0] <- NA
            cbind(table, `Pr(>Chi)` = pchisq(vals, abs(dfs), lower.tail = FALSE))
        })
    }
}

#*********************************************************************#
#**************** score.jmdem (score test for JMDEM) *****************#
#*********************************************************************#
score.jmdem <- function(object, ...)
{
    ### Determine row name by using the model names for comaprison
    tempnam <- if (class(object) == "list") {unlist(lapply(strsplit(gsub("list|\\(|\\)", "", deparse(substitute(object))), ","), trimws))} else {deparse(substitute(object))}
    
    ### Evaluate the dots
    dotargs <- list(...)
    if (length(dotargs))
    {
        if (!class(object) == "list")
        {
            object <- list(object)
        }
        object <- c(object, dotargs)
        object <- object[unlist(lapply(object, function(x) inherits(x, "jmdem")))]
        tempnam <- c(tempnam, sapply(substitute(list(...))[-1], deparse))
    }
    
    ### Determine the number of models involved
    nmodels <- length(object)
    if (nmodels == 1) {return("Test cannot be conducted because only 1 model has been specified.")}
    
    ### Execute score test
    test.res <- matrix(0, 0, 1)
    rownames.vec <- character(0)
    for (i in 1:(nmodels - 1))
    {
        for (j in (i + 1):nmodels)
        {
            ### Assign the "smaller" model to object 1
            object1 <- if (length(object[[i]]$coefficients) < length(object[[j]]$coefficients)) {object[[i]]} else {object[[j]]}
            object2 <- if (length(object[[i]]$coefficients) < length(object[[j]]$coefficients)) {object[[j]]} else {object[[i]]}
            
            ### Determine the parameter to be tested
            p0 <- vwhich(names(object1$beta), names(object2$beta), reverse = TRUE, length.check = TRUE)
            q0 <- vwhich(names(object1$lambda), names(object2$lambda), reverse = TRUE, length.check = TRUE)
            
            ### Extend the parameter vector from the "smaller" model with the tested parameters
            h0.value <- 0
            testxnum <- abs(length(object2$beta) - length(object1$beta))
            testznum <- abs(length(object2$lambda) - length(object1$lambda))
            m0 <- c(object1$beta, rep(h0.value, testxnum))
            d0 <- c(object1$lambda, rep(h0.value, testznum))
            
            ### Rank of the fitted model
            mp <- object2$mean.rank
            dp <- object2$dispersion.rank
            
            ### Create the length of the valid parameters based on the rank
            mp1 <- 1L:mp
            dp1 <- 1L:dp
            p1 <- 1L:(mp + dp)
            
            ### Create the valid indices
            mpivot <- if (mp) {object2$mean.qr$pivot[mp1]} else {numeric(0)}
            dpivot <- if (dp) {object2$dispersion.qr$pivot[dp1]} else {numeric(0)}
            model.pivot <- c(mpivot, dpivot + length(object2$beta))
            
            ### Indicies of the parameters
            beta.index <- if (length(object2$beta)) {1L:length(object2$beta)} else {numeric(0)}
            lambda.index <- if (length(object2$lambda)) {(length(object2$beta) + 1):length(object2$coefficients)} else {numeric(0)}
            xpivot <- object2$x
            zpivot <- object2$z
            
            ### Compute the score function and information matrix
            score.list <- list(par = c(m0, d0), y = object1$y, x = xpivot, z = zpivot, weights = object1$prior.weights, 
                               beta.index = beta.index, lambda.index = lambda.index, mfamily = object1$mean.family, 
                               dfamily = object1$dispersion.family, moffset = object1$mean.offset, doffset = object1$dispersion.offset, 
                               dev.type = object1$deviance.type, disp.adj = object1$dispersion.adjustment, df.adj = object1$df.adjustment, 
                               null.approx = object1$control$null.approx)
            add.list <- list(density = FALSE, full.loglik = FALSE, allobs = TRUE)
            g0 <- do.call("gloglik", append(score.list, add.list))
            i0 <- if (object1$information.type == "Fisher") {do.call("fisherinfo", score.list)} else if (object1$information.type == "Hessian") {-do.call("optimHess", c(list(fn = quote(loglik), gr = if (object1$gradient) {quote(gloglik)} else {NULL}), score.list, add.list))}
            c0 <- minv(i0)$inv
            
            ### Calculate the test statistic
            sctest <- t(g0) %*% c0 %*% g0
            test.res <- rbind(test.res, sctest)
            rownames.vec <- c(rownames.vec, paste(tempnam[i], tempnam[j], sep = " vs. "))
        }
    }
    
    colnames(test.res) <- "Rao"
    rownames(test.res) <- rownames.vec
    return(test.res)
}

#*********************************************************************#
#***************** wald.jmdem (Wald test for JMDEM) ******************#
#*********************************************************************#
wald.jmdem <- function(object, ...)
{
    ### Determine row name by using the model names for comaprison
    tempnam <- if (class(object) == "list") {unlist(lapply(strsplit(gsub("list|\\(|\\)", "", deparse(substitute(object))), ","), trimws))} else {deparse(substitute(object))}
    
    ### Evaluate the dots
    dotargs <- list(...)
    if (length(dotargs))
    {
        if (!class(object) == "list")
        {
            object <- list(object)
        }
        object <- c(object, dotargs)
        object <- object[unlist(lapply(object, function(x) inherits(x, "jmdem")))]
        tempnam <- c(tempnam, sapply(substitute(list(...))[-1], deparse))
    }
    
    ### Determine the number of models involved
    nmodels <- length(object)
    if (nmodels == 1) {return("Test cannot be conducted because only 1 model has been specified.")}
    
    ### Execute Wald test
    test.res <- matrix(0, 0, 1)
    rownames.vec <- character(0)
    for (i in 1:(nmodels - 1))
    {
        for (j in (i + 1):nmodels)
        {
            ### Assign the "smaller" model to object 1
            object1 <- if (length(object[[i]]$coefficients) < length(object[[j]]$coefficients)) {object[[i]]} else {object[[j]]}
            object2 <- if (length(object[[i]]$coefficients) < length(object[[j]]$coefficients)) {object[[j]]} else {object[[i]]}
            
            ### Determine the parameter to be tested
            h0.value <- 0
            p0 <- 1L - as.numeric(as.logical(match(names(object2$beta), names(object1$beta), nomatch = 0L)))
            q0 <- 1L - as.numeric(as.logical(match(names(object2$lambda), names(object1$lambda), nomatch = 0L)))
            par0 <- diag(c(p0, q0))
            hyp0 <- rep(h0.value, length(object2$coefficients))
            
            ### Compute the sandwich estimator of the Wald statistic
            g0 <- par0 %*% object2$coefficients - hyp0
            i0 <- minv(t(par0) %*% minv(object2$info.matrix)$inv %*% par0)$inv
            
            ### Calculate the test statistic
            waldtest <- t(g0) %*% i0 %*% g0
            test.res <- rbind(test.res, waldtest)
            rownames.vec <- c(rownames.vec, paste(tempnam[i], tempnam[j], sep = " vs. "))
        }
    }
    
    colnames(test.res) <- "Wald"
    rownames(test.res) <- rownames.vec
    return(test.res)
}

#*********************************************************************#
#****** loglik (log-likelihood function of exponential family) *******#
#*********************************************************************#
loglik <- function(par, y, x, z, weights, mfamily, dfamily, beta.index, lambda.index,
                   moffset = NULL, doffset = NULL, dev.type = c("deviance", "pearson"), 
                   disp.adj = FALSE, df.adj = FALSE, density = FALSE, full.loglik = FALSE, 
                   allobs = FALSE, null.approx = 1e-08)
{
    ### Check whether some character inputs are correct
    dev.type <- match.arg(dev.type)
    
    ### Sample size
    n <- nrow(as.matrix(y))

    ### Number of parameters
    npar <- length(par)
    
    ### Unify the specification of the deviance type
    dev.type <- tolower(dev.type)

    ### Extract mean and dispersion effects
    beta <- par[beta.index]
    lambda <- par[lambda.index]
    
    ### Check and extend call list by offset
    if (is.null(moffset)) {moffset <- rep(0L, n)}
    if (is.null(doffset)) {doffset <- rep(0L, n)}
    
    ### Compute mean and dispersion parameters
    eta <- if (ncol(x)) {x %*% beta + moffset} else {rep(0L, n)}
    delta <- if (ncol(z)) {z %*% lambda + doffset} else {rep(0L, n)}
    mu <- mfamily$linkinv(eta)
    phi <- dfamily$linkinv(delta)
    comp <- loglik.comp(family = mfamily, y = y, mu = mu, weights = weights, phi = phi)

    ### Degrees of freedom adjustment
    dfadj <- if (df.adj){(n - npar) / n} else {1}

    ### Dispersion adjustment
    if (disp.adj && dev.type == "deviance") {phi <- phi * (1 + comp$dispadj)}

    ### Core term and variance function of the likelihelood function ###
    vf <- if (dev.type == "deviance") {mfamily$variance(y)} else if (dev.type == "pearson") {mfamily$variance(mu)}
    dev <- if (dev.type == "deviance") {mfamily$dev.resids(y = y, mu = mu, wt = weights)} else if (dev.type == "pearson") {weights * (y - mu) ^ 2 / mfamily$variance(mu)}

    ### Approximation for zero-observations ###
    vf[vf == 0] <- null.approx
    
    ### Compute likelihood ###
    res <- if (!full.loglik) {-0.5 * (dfadj * log(2 * pi * phi / weights * vf) + dev / phi)} else {sum(comp$cfunc + weights / phi * (y * comp$theta - comp$btheta))}
    if (!allobs) {res <- sum(res)}
    
    ### Convert to density if required ###
    if (density) {res <- exp(res)}

    return(res)
}

#*********************************************************************#
#************************* gradient loglik  **************************#
#*********************************************************************#
gloglik <- function(par, y, x, z, weights, mfamily, dfamily, beta.index, lambda.index, 
                    moffset = NULL, doffset = NULL, dev.type = c("deviance", "pearson"), 
                    disp.adj = FALSE, df.adj = FALSE, density = FALSE, full.loglik = FALSE, 
                    allobs = FALSE, null.approx = 1e-08)
{
    ### Check whether some character inputs are correct
    dev.type <- match.arg(dev.type)

    ### Sample size    
    n <- nrow(as.matrix(y))
    
    ### Number of parameters
    npar <- length(par)
    
    ### Unify the specification of the deviance type
    dev.type <- tolower(dev.type)
    
    ### Extract mean and dispersion effects
    beta <- par[beta.index]
    lambda <- par[lambda.index]

    ### Check and extend call list by offset
    if (is.null(moffset)) {moffset <- rep(0L, n)}
    if (is.null(doffset)) {doffset <- rep(0L, n)}
    
    ### Compute mean and dispersion parameters
    eta <- if (ncol(x)) {x %*% beta + moffset} else {rep(0L, n)}
    delta <- if (ncol(z)) {z %*% lambda + doffset} else {rep(0L, n)}
    mu <- mfamily$linkinv(eta)
    phi <- dfamily$linkinv(delta)
    comp <- loglik.comp(family = mfamily, y = y, mu = mu, weights = weights, phi = phi)

    ### Compute the first and second cumulant orders of the parameters
    vmu <- mfamily$variance(mu)
    gmu <- mfamily$mu.eta(eta)
    vphi <- dfamily$variance(phi)
    gphi <- dfamily$mu.eta(delta)
    vmudmu <- comp$vmudmu
    
    ### Check if the model is nested
    mnest <- length(which(names(beta) %in% c("delta"))) > 0
    dnest <- length(which(names(lambda) %in% c("eta"))) > 0

    ### Degrees of freedom adjustment
    dfadj <- if (df.adj){(n - npar) / n} else {1}
    
    ### Dispersion adjustment
    if (disp.adj && dev.type == "deviance") 
    {
        phi <- phi * (1 + comp$dispadj)
        vphi <- vphi * dfamily$variance(1 + comp$dispadj)
    }

    ### Calculate deviance
    dev <- if (dev.type == "deviance") {mfamily$dev.resids(y = y, mu = mu, wt = weights)} else if (dev.type == "pearson") {dev <- weights * (y - mu) ^ 2 / vmu}
    
    ### Derivatives of the constraint parts
    psix <- if (dnest) {lambda[which(names(lambda) %in% "eta")] * x} else {as.matrix(0 * x)}
    psiz <- if (mnest) {beta[which(names(beta) %in% "delta")] * z} else {as.matrix(0 * z)}
    
    ### Construct score function parts
    grfx0 <- t(x) %*% diag(as.vector((weights * gmu) / (phi * vmu))) %*% (y - mu)
    grfx1 <- if (dev.type == "pearson") {t(x) %*% diag(as.vector((vmudmu * gmu) / (2 * phi * vmu))) %*% (dev - phi)} else {as.matrix(0 * beta)}
    grfx2 <- t(psix) %*% diag(as.vector(gphi / (2 * vphi))) %*% (dev - phi)
    grfz0 <- t(z) %*% diag(as.vector(gphi / (2 * vphi))) %*% (dev - phi)
    grfz1 <- t(psiz) %*% diag(as.vector((weights * gmu) / (phi * vmu))) %*% (y - mu)
    
    ### Add score function together
    grfx <- grfx0 + grfx1 + grfx2
    grfz <- grfz0 + grfz1
    
    return(c(grfx, grfz))
}

#*********************************************************************#
#******************** Fisher information matrix **********************#
#*********************************************************************#
fisherinfo <- function(par, x, y, z, weights, moffset, doffset, beta.index, lambda.index, mfamily, 
                       dfamily, dev.type = c("deviance", "pearson"), disp.adj = FALSE, df.adj = FALSE, 
                       null.approx = 1e-08)
{
    ### Check whether some character inputs are correct
    dev.type <- match.arg(dev.type)
    
    ### Sample size    
    n <- nrow(as.matrix(y))
    
    ### Number of parameters
    npar <- length(par)
    
    ### Unify the specification of the deviance type
    dev.type <- tolower(dev.type)
    
    ### Extract mean and dispersion effects
    beta <- par[beta.index]
    lambda <- par[lambda.index]
    
    ### Check and extend call list by offset
    if (is.null(moffset)) {moffset <- rep(0L, n)}
    if (is.null(doffset)) {doffset <- rep(0L, n)}
    
    ### Compute mean and dispersion parameters
    eta <- x %*% beta + moffset
    delta <- z %*% lambda + doffset
    mu <- mfamily$linkinv(eta)
    phi <- dfamily$linkinv(delta)
    comp <- loglik.comp(family = mfamily, y = y, mu = mu, weights = weights, phi = phi)
    
    ### Compute the first and second cumulant orders of the parameters
    vmu <- mfamily$variance(mu)
    gmu <- mfamily$mu.eta(eta)
    vphi <- dfamily$variance(phi)
    gphi <- dfamily$mu.eta(delta)
    vmudmu <- comp$vmudmu
    
    ### Check if the model is nested
    mnest <- length(which(names(beta) %in% c("delta"))) > 0
    dnest <- length(which(names(lambda) %in% c("eta"))) > 0

    ### Degrees of freedom adjustment
    dfadj <- if (df.adj){(n - npar) / n} else {1}

    ### Dispersion adjustment
    if (disp.adj && dev.type == "deviance") 
    {
        phi <- phi * (1 + comp$dispadj)
        vphi <- vphi * dfamily$variance(1 + comp$dispadj)
    } 

    ###  Calculate deviance
    dev <- if (dev.type == "deviance") {mfamily$dev.resids(y = y, mu = mu, wt = weights)} else if (dev.type == "pearson") {dev <- weights * (y - mu) ^ 2 / vmu}
    
    ### Functions for reweighted matrix
    w.m <- diag(as.vector((weights * gmu ^ 2) / (vmu * phi)))
    v.m <- diag(as.vector(gphi ^ 2/ (2 * vphi)))
    t.m <- diag(as.vector((vmudmu * gmu) ^ 2 / (2 * vmu ^ 2)))
    u.m <- diag(as.vector((vmudmu * gmu * gphi) / (2 * phi * vmu)))
    z.m <- diag(as.vector(0 * mu))
    
    ### Constraints
    psix <- if (dnest) {lambda[which(names(lambda) %in% "eta")] * x} else {as.matrix(0 * x)}
    psiz <- if (mnest) {beta[which(names(beta) %in% "delta")] * z} else {as.matrix(0 * z)}

    ### Blocks of expected Hessian (Used for both effects and deviance types)
    xwx <- t(x) %*% w.m %*% x
    zvz <- t(z) %*% v.m %*% z
    xwpz <- t(x) %*% w.m %*% psiz
    zvpx <- t(z) %*% v.m %*% psix
    pzwpz <- t(psiz) %*% w.m %*% psiz
    pxvpx <- t(psix) %*% v.m %*% psix
    block0 <- t(x) %*% z.m %*% z
    if (dev.type == "deviance") 
    {
        # Block elements
        ibb <- if (dnest) {xwx + pxvpx} else {xwx} # Top left
        ibg <- if (mnest) {xwpz} else if (dnest) {t(zvpx)} else {block0} # Top right and bottom left
        igg <- if (mnest) {zvz + pzwpz} else {zvz} # Bottom right
    } 
    else if (dev.type == "pearson") 
    {
        # Rest of the block components
        xtx <- t(x) %*% t.m %*% x
        xux <- t(x) %*% u.m %*% x
        xuz <- t(x) %*% u.m %*% z
        xwpx <- t(x) %*% w.m %*% psix
        xtpz <- t(x) %*% t.m %*% psiz
        ztpz <- t(z) %*% t.m %*% psiz
        pztpz <- t(psiz) %*% t.m %*% psiz
        x2ux <- t(x) %*% (2 * u.m) %*% x
        z2upz <- t(z) %*% (2 * u.m) %*% psiz
        # Block elements
        ibb <- if (dnest) {xwx + xtx + x2ux + pxvpx} else {xwx + xtx} # Top left
        ibg <- if (mnest) {xuz + xwpz + xtpz} else if (dnest) {xuz + t(zvpx)} else {xuz} # Top right and bottom left
        igg <- if (mnest) {zvz + z2upz + pzwpz + pztpz} else {zvz} # Bottom right
    }
    return(rbind(cbind(ibb, ibg), cbind(t(ibg), igg)))
}

#*********************************************************************#
#********** convert estimated mean to canonical parameter ************#
#*********************************************************************#
loglik.comp <- function(family, y, mu, weights, phi)
{
    ### Sample size
    n <- nrow(as.matrix(y))
    w <- weights / phi

    ### Compute the components based on the distribution
    if (family$family == "gaussian")
    {
        theta <- mu
        btheta <- 1 / 2 * (theta ^ 2)
        btheta1 <- theta
        cfunc <- -0.5 * (y ^ 2 * w + log(2 * pi / w))
        dispadj <- 0 #b(phi, mu)
        cumadj <- rep(1, n) #1+rho4/2
        vmudmu <- rep(0, n)
    }
    else if (family$family == "inverse.gaussian")
    {
        theta <- -1 / (2 * (mu ^ 2))
        btheta <- -((-2 * theta) ^ 0.5)
        btheta1 <- (-2 * theta) ^ (-0.5)
        cfunc <- -0.5 * (log(2 * pi / w * y ^ 3) + (w / y))
        dispadj <- 0
        cumadj <- 1 + 15 * phi * mu / 2
        vmudmu <- 3 * mu ^ 2
    }
    else if (family$family == "Gamma")
    {
        theta <- -1 / mu
        btheta <- -log(-theta)
        btheta1 <- -1 / theta
        cfunc <- w * log(w * y) - log(y) - lgamma(w)
        dispadj <- phi / 6
        cumadj <- 1 + 3 * phi
        vmudmu <- 2 * mu
    }
    else if (family$family == "binomial")
    {
        theta <- log(mu / (1 - mu))
        btheta <- log(1 + exp(theta))
        btheta1 <- exp(theta) / (1 + exp(theta))
        cfunc <- log(w)
        dispadj <- (phi / (6 * weights)) * ((1 - mu * (1 - mu)) / (mu * (1 - mu)))
        cumadj <- 1 + (phi / (2 * weights)) * ((1 - 6 * mu * (1 - mu)) / (mu * (1 - mu)))
        vmudmu <- (1 - 2 * mu) / weights
    }
    else if (family$family == "poisson")
    {
        theta <- log(mu)
        btheta <- exp(theta)
        btheta1 <- exp(theta)
        cfunc <- -lfactorial(y)
        dispadj <- phi / (6 * mu)
        cumadj <- 1 + phi / (2 * mu)
        vmudmu <- rep(1, n)
    }
    return(list(theta = theta, btheta = btheta, btheta1 = btheta1, cfunc = cfunc, 
                dispadj = dispadj, cumadj = cumadj, vmudmu = vmudmu))
}

#*********************************************************************#
#******************* jmdem.sim (JMDEM simulations) *******************#
#*********************************************************************#
jmdem.sim <- function(mformula = "y ~ 1 + x", dformula = "~ 1 + z", data = NULL, 
                      beta.true, lambda.true, mfamily = gaussian, 
                      dfamily = Gamma, dev.type = c("deviance", "pearson"), 
                      x.str = list(type = "numeric", random.func = "runif", param = list()), 
                      z.str = list(type = "numeric", random.func = "runif", param = list()), 
                      n = NULL, simnum = NULL, trace = FALSE, asymp.test = FALSE, 
                      weights = NULL, moffset = NULL, doffset = NULL, 
                      mustart = NULL, phistart = NULL, betastart = NULL, 
                      lambdastart = NULL, hessian = TRUE, na.action, 
                      grad.func = TRUE, fit.method = "jmdem.fit", 
                      method = c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN", "Brent"), 
                      df.adj = FALSE, disp.adj = FALSE, full.loglik = FALSE, 
                      mcontrasts = NULL, dcontrasts = NULL, beta.first = TRUE, 
                      prefit = TRUE, control = list(...), 
                      minv.method = c("solve", "chol2inv", "ginv"), ...)
{
    
    ### Listing of all model arguments
    orgmatch <- match.call()
    arglist <- as.list(orgmatch)[-1]
    argnam <- names(arglist)
    simcall <- orgmatch
    
    ### Check whether some character inputs are correct
    dev.type <- match.arg(dev.type)
    method <- match.arg(method)
    minv.method <- match.arg(minv.method, several.ok = TRUE)
    
    ### Check whether intercepts are included in the models
    mintercept <- as.logical(attr(terms(mformula), "intercept"))
    dintercept <- as.logical(attr(terms(dformula), "intercept"))
    
    ### Create the function argument to create data
    arg.jmdem.sim.crdat <- as.list(args(simdata.jmdem.sim))
    arg.jmdem.sim.crdat <- arg.jmdem.sim.crdat[!names(arg.jmdem.sim.crdat) == "" & !names(arg.jmdem.sim.crdat) == "..."]
    crdata.args <- lwhich(arglist, arg.jmdem.sim.crdat)
    
    ### Calculate number of total parameters
    xnum <- as.numeric(mintercept) + length(attr(terms(mformula), "term.labels"))
    znum <- as.numeric(dintercept) + length(attr(terms(dformula), "term.labels"))

    ### Sample size
    if (is.null(n)) 
    {
        n <- if (is.null(data)) {return("Sample size not specified!")} else {nrow(data[[1]])}
    }

    ### Simulation runs
    if (is.null(simnum)) 
    {
        simnum <- if (is.null(data)) {return("Number of simulation runs not specified!")} else {length(data)}
    }
    
    ### Observation weights
    if (is.null(weights)) {weights <- rep(1, n)}
    
    ### Start simulations
    simlist <- list()
    for (j in 1:simnum)
    {
        # Create data if no user data are provided
        if (is.null(data)) 
        {
            simdata <- do.call("simdata.jmdem.sim", crdata.args)
            arg.sim <- if (simnum == 1) {list(data = simdata)} else {list(data = simdata[[1]])}
        }    
        else
        {
            simdata <- eval(data[[j]])
            arglist$data <- simdata
            arg.sim <- list()
        }
        
        ### Prepare user arguments for fitting
        arg.jmdem <- as.list(args(jmdem))
        arg.sim <- c(lwhich(arglist, arg.jmdem), arg.sim)
        
        ### Fit a model by simulation data
        simres <- try(do.call("jmdem", arg.sim), silent = TRUE)
        
        ### Complete this simulation only if simulation was successful
        if (!class(simres) == "try-error") 
        {
            names(simres$coefficients) <- c(paste("beta", names(simres$coefficients)[1:xnum], sep = "."), paste("lambda", names(simres$coefficients)[-(1:xnum)], sep = "."))
            names(simres$coefficients)[grep("(Intercept)", names(simres$coefficients))] <- c(if (mintercept) {"beta.x0"}, if (dintercept) {"lambda.z0"})
            
            ### Carry out asymptotic tests
            if (asymp.test)
            {
                ### Parameter for asynptotic tests
                df.test <- simres$df.null - simres$df.residual # Degrees of freedom
                arg.sim0 <- arg.sim
                arg.sim0$mformula <- as.formula("y ~ 1")
                arg.sim0$dformula <- as.formula(" ~ 1")
                simnull <- do.call("jmdem", arg.sim0) # Null model
                
                ### Asynptotic tests
                waldx <- wald.jmdem(simres, simnull) # Wald test
                scorex <- score.jmdem(simres, simnull) # Score test
            }
            
            ### Return simulation results
            curr.sim <- append(simres, list(simcall = simcall, beta.true = beta.true, lambda.true = lambda.true, asymp.test = asymp.test))
            if (asymp.test) {curr.sim <- append(curr.sim, list(wald = waldx, rao = scorex))}
            simlist <- append(simlist, list(curr.sim))
            
            ### Print current iteration estimates if trace is TRUE
            if (trace)
            {
                trace.out <- t(as.matrix(c(simres$coefficients, if (asymp.test) {c(wald = waldx, rao = scorex)} else {NULL})))
                rownames(trace.out) <- j
                print(trace.out)
            }
        }
            
    }
    
    class(simlist) <- c("jmdem.sim")
    return(simlist)
}

#*********************************************************************#
#****** simdata.jmdem.sim (generate data for JMDEM simulations) ******#
#*********************************************************************#
simdata.jmdem.sim <- function(mformula = "y ~ 1 + x", dformula = "~ 1 + z", beta.true, lambda.true, 
                              x.str = list(type = "numeric", random.func = "runif", param = list()), 
                              z.str = list(type = "numeric", random.func = "runif", param = list()), 
                              mfamily = gaussian, dfamily = Gamma, weights = NULL, n, simnum = 1, 
                              moffset = NULL, doffset = NULL)
{
    
    ### Listing of all model arguments
    orgmatch <- match.call()
    arglist <- as.list(orgmatch)[-1]
    argnam <- names(arglist)
    simcall <- orgmatch
    
    ### Convert formula
    mformula <- as.formula(mformula)
    dformula <- as.formula(dformula)
    
    ### Convert distribution family
    mfamily <- check.family(family = mfamily)
    dfamily <- check.family(family = dfamily, type = "dispersion")
    
    ### Fill in default setting options
    if (is.null(arglist$mformula)) {arglist$mformula <- mformula}
    if (is.null(arglist$dformula)) {arglist$dformula <- dformula}
    if (is.null(arglist$mfamily)) {arglist$mfamily <- mfamily}
    if (is.null(arglist$dfamily)) {arglist$dfamily <- dfamily}
    
    ### Variable name of Y
    ynames <- as.character(attr(terms(mformula), "variables")[[2]])
    
    ### Variable name of X and Z
    xnames <- attr(terms(mformula), "term.labels")
    znames <- attr(terms(dformula), "term.labels")
    
    ### Intercept term
    xicept <- attr(terms(mformula), "intercept")
    zicept <- attr(terms(dformula), "intercept")
    
    ### Vector of effect orders
    xorder <- attr(terms(mformula), "order")
    zorder <- attr(terms(dformula), "order")

    ### Calculate number of total parameters
    xnum <- xicept + length(xnames)
    znum <- zicept + length(znames)
    
    ### Create data 
    # Check and adjust user settings for the design matrices
    strlist.default <- list(type = "numeric", random.func = "runif", param = list())
    x.misslist <- lwhich(strlist.default, x.str, reverse = TRUE)
    z.misslist <- lwhich(strlist.default, z.str, reverse = TRUE)
    if (length(x.misslist)) {x.str <- c(x.str, x.misslist)}
    if (length(z.misslist)) {z.str <- c(z.str, z.misslist)}
    
    ### Prepare conversion function for the predictors
    xtype <- paste("as", rep(x.str$type, length(xnames)), sep = ".")
    ztype <- paste("as", rep(z.str$type, length(znames)), sep = ".")
    
    ### Prepare random function for the predictors
    xfun <- rep(x.str$random.func, length(xnames))
    zfun <- rep(z.str$random.func, length(znames))
    
    ### Observation weights
    if (is.null(weights)) {weights <- rep(1, n)}
    
    ### Output list
    datlist <- list()
    
    ### Generate Data
    for (j in 1:simnum)
    {
        # Design Matrix X
        tempdat <- data.frame(matrix(0, n, 0))
        for (i in 1:length(xnames))
        {
            if (xorder[i] == 1)
            {
                tempvar <- do.call(xfun[i], args = c(list(n = n), x.str$param))
                if (xtype[i] == "as.factor") {tempvar <- round(tempvar, 0)}
                tempvar <- as.data.frame(do.call(xtype[i], list(x = tempvar)))
                tempdat <- cbind(tempdat, tempvar)
            }        
        }
        colnames(tempdat) <- xnames[xorder == 1]
        simdata <- tempdat
        
        # Design Matrix Z
        tempdat <- data.frame(matrix(0, n, 0))
        for (i in 1:length(znames))
        {
            if (zorder[i] == 1)
            {
                tempvar <- do.call(zfun[i], args = c(list(n = n), z.str$param))
                if (ztype[i] == "as.factor") {tempvar <- round(tempvar, 0)}
                tempvar <- as.data.frame(do.call(ztype[i], list(x = tempvar)))
                tempdat <- cbind(tempdat, tempvar)
            }
        }
        colnames(tempdat) <- znames
        simdata <- cbind(simdata, tempdat)

        # Set offsets
        if (is.null(moffset)) {moffset <- rep(0, n)}
        if (is.null(doffset)) {doffset <- rep(0, n)}
        
        # Add delta or eta if models are nested 
        if (match("delta", xnames, 0L)) 
        {
            temp.z <- as.matrix(cbind(matrix(1, n, zicept), simdata[, znames[!znames == "eta"], drop = FALSE]))
            znames.temp <- c(if (zicept) {"(intercept)"}, znames)
            simdata$delta <- temp.z %*% as.matrix(lambda.true[!znames.temp == "eta"]) + doffset
        }
        if (match("eta", znames, 0L)) 
        {
            temp.x <- as.matrix(cbind(matrix(1, n, xicept), simdata[, xnames[!xnames == "delta"], drop = FALSE]))
            xnames.temp <- c(if (xicept) {"(intercept)"}, xnames)
            simdata$eta <- temp.x %*% as.matrix(beta.true[!xnames.temp == "delta"]) + moffset
        }

        # Response Y
        y <- data.frame(rep(0, n))
        colnames(y) <- ynames
        simdata <- cbind(simdata, y)
        
        # Compute individual response mean and dispersion
        x <- model.matrix(mformula, simdata)
        z <- model.matrix(dformula, simdata)
        mu <- mfamily$linkinv(x %*% beta.true + moffset)
        phi <- dfamily$linkinv(z %*% lambda.true + doffset)
        variance <- mfamily$variance(mu) * phi
        
        # Generate response variable following the underlying distribution and parameters 
        if (mfamily$family == "gaussian")
        {
            y <- rnorm(n = n, mean = mu, sd = variance ^ 0.5)
        }
        else if (mfamily$family == "inverse.gaussian")
        {
            requireNamespace("statmod", quietly = TRUE)
            y <- statmod::rinvgauss(n = n, mean = mu, dispersion = variance / (mu ^ 3))
        }
        else if (mfamily$family == "Gamma")
        {
            y <- rgamma(n = n, scale = variance / mu, shape = mu ^ 2 / variance)
        }
        else if (mfamily$family == "binomial")
        {
            requireNamespace("VGAM", quietly = TRUE)
            new.weights <- weights
            new.weights[weights - 1 == 0] <- 1.1
            rho <- variance / (new.weights * mu * (1 - mu) * (new.weights - 1)) - 1 / (new.weights - 1)
            y <- VGAM::rbetabinom(n = n, size = weights, prob = mu, rho = rho)
        }
        else if (mfamily$family == "poisson")
        {
            k <- mu / (phi - 1)
            y <- rnbinom(n = n, size = k, mu = mu)
        }
        simdata[, ynames] <- y
        datlist <- append(datlist, list(simdata))
    }
    
    return(if (simnum == 1) {datlist[[1]]} else {datlist})
}

#*********************************************************************#
#********** getdata.jmdem.sim (save JMDEM simulation data) ***********#
#*********************************************************************#
getdata.jmdem.sim <- function(object)
{
    if (!class(object) == "jmdem.sim") {return("The function is only designed for jmdem.sim objects.")}
    
    simdat <- list()
    for(i in 1:length(object))
    {
        simdat <- append(simdat, list(object[[i]]$data))
    }
    
    return(simdat)
}

#*********************************************************************#
#*********** summary.jmdem.sim (JMDEM simulation summary) ************#
#*********************************************************************#
summary.jmdem.sim <- function(object, digits = max(3L, getOption("digits") - 3L), 
                              scientific = FALSE, pvalue = 0.05, 
                              minv.method = c("solve", "chol2inv", "ginv"), 
                              other.call = FALSE, details = FALSE, ...)
{
    ### Listing of all function arguments
    orgmatch <- match.call()
    arglist <- as.list(orgmatch)[-1]
    
    ### Check whether some character inputs are correct
    minv.method <- match.arg(minv.method, several.ok = TRUE)
    
    ### Set default digits and scientific notation settings
    digits <- if (!is.null(arglist$digits)) {digits} else {max(3L, getOption("digits") - 3L)}
    scientific <- if (!is.null(arglist$scientific)) {scientific} else {FALSE}

    ### Row and column names
    rownum <- if (length(object) > 1) {1:length(object)} else {1}
    
    ### Coefficient summary
    coefsum <- t(as.data.frame(lapply(object, coef)))
    rownames(coefsum) <- rownum
    colnames(coefsum) <- paste("coef", names(coef(object[[1]])), sep = ".")

    ### Standard error summary
    sesum <- t(as.data.frame(lapply(object, std.err, include.covmat = FALSE)))
    rownames(sesum) <- rownum
    colnames(sesum) <- paste("stderr", names(coef(object[[1]])), sep = ".")

    ### True value summary
    betatruesum <- lapply(object, '[[', "beta.true")
    lambdatruesum <- lapply(object, '[[', "lambda.true")
    truesum <- t(as.data.frame(mapply(c, betatruesum, lambdatruesum, SIMPLIFY=FALSE)))
    
    ### Degrees of freedom summary
    dfsum <- t(as.data.frame(lapply(object, '[[', "df.residual")))
    
    ### Confidence interval coverage
    cilowersum <- coefsum + qt(rep(pvalue/2, nrow(coefsum)), df = dfsum, lower.tail = TRUE) * sesum
    ciuppersum <- coefsum + qt(rep(pvalue/2, nrow(coefsum)), df = dfsum, lower.tail = FALSE) * sesum
    colnames(ciuppersum) <- names(coef(object[[1]]))
    colnames(cilowersum) <- names(coef(object[[1]]))
    coverage <- as.data.frame(truesum >= cilowersum & truesum <= ciuppersum)
    rownames(coverage) <- rownum
    colnames(coverage) <- paste("cover", names(coef(object[[1]])), sep = ".")
    
    ### Confidence interval summary
    cisumcolnum <- order(c(seq(1, ncol(cilowersum) * 2 - 1, by = 2), seq(2, ncol(cilowersum) * 2, by = 2)))
    cisum <- as.matrix(cbind(cilowersum, ciuppersum)[, cisumcolnum])
    if (length(object) == 1) {cisum <- t(cisum)}
    rownames(cisum) <- rownum
    colnames(cisum) <- paste(c("cil", "ciu"), colnames(cisum), sep = ".")
    
    ### Iteration number summary
    itersum <- t(as.data.frame(lapply(object, '[[', "iter")))
    colnames(itersum) <- "iterations"
    rownames(itersum) <- rownum
    colnames(itersum) <- "iterations"

    ### Asymptotic test summary
    atest <- object[[1]]$asymp.test
    if (atest)
    {
        asymptsum <- cbind(t(as.data.frame(lapply(object, '[[', "wald"))), t(as.data.frame(lapply(object, '[[', "rao"))))
        rownames(asymptsum) <- rownum
        colnames(asymptsum) <- c(" Wald-Test", " Rao-Test")
        avasympttab <- t(as.data.frame(colMeans(asymptsum)))
        rownames(avasympttab) <- c("Chi-square values ")
    }
    
    ### Summarise results
    avrestab <- rbind(colMeans(coefsum), colMeans(sesum), colMeans(coverage))
    colnames(avrestab) <- gsub("coef.", " ", colnames(coefsum))
    rownames(avrestab) <- c("coefficients", "std. errors", "CI coverage")
    
    ans <- list(digits = digits, scientific = scientific, details = details, other.call = other.call, pvalue = pvalue, 
                beta.true = betatruesum[[1]], lambda.true = lambdatruesum[[1]], simcall = object[[1]]$simcall, 
                mformula = object[[1]]$mean.formula, dformula = object[[1]]$dispersion.formula, 
                mfamily = object[[1]]$mean.family, dfamily = object[[1]]$dispersion.family, 
                coefficients = coefsum, stderr = sesum, iterations = itersum, confint = cisum, 
                coverage = coverage, average.summary = avrestab, asymp.test = atest)
    if (atest) {ans <- c(ans, list(wald.rao.test = asymptsum, average.asymp.test = avasympttab))}
    
    class(ans) <- "summary.jmdem.sim"
    return(ans)
}

#*********************************************************************#
#****** print.summary.jmdem.sim (print JMDEM simulation summary) *****#
#*********************************************************************#
print.summary.jmdem.sim <- function (x, digits = max(3L, getOption("digits") - 3L), scientific = FALSE, 
                                     pvalue = 0.05, signif.stars = getOption("show.signif.stars"), 
                                     other.call = FALSE, details = FALSE, ...)
{
    ### Edit mean formula
    mform.char <- as.character(x$mformula)
    mform.terms <- terms(x$mformula)
    mintercept <- as.logical(attr(mform.terms, "intercept"))
    mform.labels <- attr(mform.terms,"term.labels")
    mform.micept <- if (mintercept) {x$beta.true[1]} else {""}
    if (mintercept) {x$beta.true <- x$beta.true[-1]}
    mform.items <- paste(paste(x$beta.true, mform.labels, sep = " * "), collapse = " + ")
    mform.edit <- paste(mform.char[2], " ", mform.char[1], " ", mform.micept, if (mintercept) {" + "}, mform.items, sep = "")

    ### Edit dispersion formula
    dform.char <- as.character(x$dformula)
    dform.terms <- terms(x$dformula)
    dintercept <- as.logical(attr(dform.terms, "intercept"))
    dform.labels <- attr(dform.terms,"term.labels")
    dform.micept <- if (dintercept) {x$lambda.true[1]} else {""}
    if (dintercept) {x$lambda.true <- x$lambda.true[-1]}
    dform.items <- paste(paste(x$lambda.true, dform.labels, sep = " * "), collapse = " + ")
    dform.edit <- paste("D ", dform.char[1], " ", dform.micept, if (dintercept) {" + "}, dform.items, sep = "")
    
    ### Listing of all function arguments
    orgmatch <- match.call()
    arglist <- as.list(orgmatch)[-1]
    
    ### Set default digits and scientific notation settings
    digits <- x$digits
    scientific <- x$scientific
    details <- x$details
    pvalue <- x$pvalue
    other.call <- x$other.call

    ### Summarise results
    restab <- cbind(x$coefficient, x$stderr)
    restab.dignum <- floor(log10(max(abs(restab)))) + digits
    restab <- format(round(restab, digits), digits = restab.dignum, scientific = scientific)
    
    ### Prepare other call options output
    if (other.call)
    {
        funnam <- x$simcall[[1]]
        x$simcall <- noquote(gsub("  ", " ", substring(gsub(".{1}$", "", gsub(funnam, "", deparse(x$simcall[grep("formula|data|family|true|simnum", names(x$simcall), invert = TRUE)], width.cutoff = 500L))), 2)))
    }
    
    ### Prepare estimation summary table
    sumtab.dignum <- floor(log10(max(x$average.summary))) + digits
    sumtab <- format(round(x$average.summary, digits), digits = sumtab.dignum, scientific = scientific, justify = "right")
    
    ### Prepare asymptotic test table
    atest <- x$asymp.test
    if (atest)
    {
        asymptest.dignum <- floor(log10(max(x$wald.rao.test))) + digits
        asymptest <- format(round(x$wald.rao.test, digits), digits = asymptest.dignum, scientific = scientific)
        
        ### Prepare asymptotic test summary table
        asymtab.dignum <- floor(log10(max(x$average.asymp.test))) + digits
        asymtab <- format(round(x$average.asymp.test, digits), digits = asymtab.dignum, scientific = scientific)
    }

    ### Output results
    cat("Summary of Simulation Studies\n")
    cat("\nMean family: ", x$mfamily$family, "(link = ", x$mfamily$link, ")", sep = "")
    cat("\nDispersion family: ", x$dfamily$family, "(link = ", x$dfamily$link, ")\n", sep = "")
    cat("\nModel (with true parameters): ")
    cat(mform.edit, "\n", sep = "") 
    cat(rep(" ", 30), dform.edit, "\n", sep = "")
    if (other.call)
    {
        cat("\nOther call options: ", "\n", x$simcall, "\n", sep = "")
    }
    cat("\nAveraged estimates of the parameters and asymptotic test statistics after", nrow(x$coefficients), "simulations:\n\n")
    print.default(sumtab, quote = FALSE, right = TRUE)
    if (atest) 
    {
        cat("\n")
        print.default(asymtab, quote = FALSE, right = TRUE)
    }
    cat("\nAvergae iterations:", mean(x$iteration))
    if (details)
    {
        cat("\n\nEstimated coefficients and the corresponding standard errors of each simulation:\n\n")
        print.default(restab, quote = FALSE, right = TRUE)
        cat("\nTrue parameter CI-Coverage of each simulation:\n\n")
        print(x$coverage, digits = digits)
        if (atest) 
        {
            cat("\nAsymptotic tests:\n\n")
            print.default(asymptest, quote = FALSE, right = TRUE)
        }
    }
}

#*********************************************************************#
#*************************** check family ****************************#
#*********************************************************************#
check.family <- function(family, type = "mean")
{
    ### Check whether some character inputs are correct
    type <- match.arg(type, choices = c("mean", "dispersion"))

        # Set default value if parameters were not specified by user 
    if (type == "dispersion")
    {
        default.dist <- stats::Gamma(link = "log")
        default.lnkfunc <- list("log")
    }
    else if (!type == "dispersion")
    {
        default.dist <- stats::gaussian()
        default.lnkfunc <- list()
    }
    
    if (is.character(family)) # If family is put inside quotation marks
    {
        if (length(grep("\\(|\\)", family)) > 0L)
        {
            famsplit <- unlist(strsplit(gsub("\"|\'|link = |link=", "", family), "\\(|\\)"))
            famfunc <- famsplit[1]
            if (famsplit[2] == ""){lnkfunc <- default.lnkfunc} else {lnkfunc <- list(famsplit[2])}
        } 
        else 
        {
            famfunc <- family
            lnkfunc <- default.lnkfunc
        }
        family <- do.call(famfunc, args = lnkfunc)
    }
    else if (is.function(family)) # If family is put as a function, i.e. no parenthesis with arguments
    {
        if (type == "dispersion"){family <- family(link = "log")} else {family <- family()}
    }
    else if (is.null(family)) # If user did not specify family at all
    {
        family <- default.dist
    }
    return(family)
}
    
#*********************************************************************#
#***************** std.err (extract standard errors) *****************#
#*********************************************************************#
std.err <- function(object, minv.method = c("solve", "chol2inv", "ginv"), range = NULL, include.covmat = TRUE, ...)
{
    ### Get the valid method to invert matrices
    if (!is.null(object$minv.method)) 
    {
        minv.method <- object$minv.method
    }
    
    #Check whether some character inputs are correct 
    if (!is.null(minv.method))
    {
        minv.method <- match.arg(minv.method)
    }
    
    ### Get the valid coefficients and its covariance matrix
    fishmat <- minv(object$info.matrix, minv.method = minv.method)
    if (is.null(range)) {range <- 1:min(dim(object$info.matrix))}
    covmat <- fishmat$inv[range, range, drop = FALSE]
    
    ### Compute the inference statistics of the parameters
    s.err <- c(list(s.err = sqrt(diag(covmat))), if (include.covmat) {list(covmat = covmat, covmat.method = fishmat$inv.nam)})
    return(s.err)
    
}

#*********************************************************************#
#****** vwhich (return vector 1 elements of found in vector 2  *******#
#*********************************************************************#
vwhich <- function(vec1, vec2, type = "vector", reverse = FALSE, length.check = FALSE)
{
    ### Unify the type argument
    type <- match.arg(arg = substr(tolower(type), 1, 1), choices = c("v", "e", "w", "b"))
    
    ### Create the element vector for the return value
    if (length.check) 
    {
        object1 <- if (length(vec1) >= length(vec2)) {vec1} else {vec2}
        object2 <- if (length(vec1) >= length(vec2)) {vec2} else {vec1}
    } 
    else 
    {
        object1 <- vec1
        object2 <- vec2
    }
    coelement <- which(vec1 %in% vec2)
    if (reverse)
    {
        vseq <- 1L:length(object1)
        coelement <- vseq[-coelement]
    }
    
    ### Create the return value based on the type
    res <- switch(type, 
                  v = object1[coelement],
                  e =, v =, w = coelement,
                  b = length(coelement) > 0,
                  "Unknown return type!")
    return(res)
}

#*********************************************************************#
#******** lwhich (return list 1 elements of found in list 2  *********#
#*********************************************************************#
lwhich <- function(list1, list2, type = "list", reverse = FALSE, length.check = FALSE)
{
    ### Unify the type argument
    type <- match.arg(arg = substr(tolower(type), 1, 1), choices = c("l", "n", "e", "v", "w", "b"))

    ### Create the element vector for the return value
    if (length.check) 
    {
        object1 <- if (length(list1) >= length(list2)) {list1} else {list2}
        object2 <- if (length(list1) >= length(list2)) {list2} else {list1}
    } 
    else 
    {
        object1 <- list1
        object2 <- list2
    }
    coelement <- which(names(object1) %in% names(object2))
    if (reverse)
    {
        lsseq <- 1L:length(object1)
        coelement <- lsseq[-coelement]
    }

    ### Create the return value based on the type
    res <- switch(type, 
                  l = object1[coelement],
                  n = names(object1)[coelement],
                  e =, v =, w = coelement,
                  b = length(coelement) > 0,
                  "Unknown return type!")
    return(res)
    
}

#*********************************************************************#
#***************** generalised inverse of a matrix  ******************#
#*********************************************************************#
ginv <- function(m, tol = sqrt(.Machine$double.eps))
{
    ### QR inverse
    dnx <- dimnames(m)
    if (is.null(dnx)) {dnx <- vector("list", 2)}
    s <- svd(m)
    nz <- s$d > tol * s$d[1]
    ans <- structure(if (any(nz)) {s$v[, nz] %*% (t(s$u[, nz])/s$d[nz])} else {m}, dimnames = dnx[2:1])
    return(ans)
}
    
#*********************************************************************#
#************************* invert a matrix  **************************#
#*********************************************************************#
minv <- function(m, minv.method = NULL, tol = sqrt(.Machine$double.eps))
{
    if (!is.null(minv.method)) 
    {
        ans <- try(do.call(minv.method, list(m)), silent = TRUE)
        if (minv.method == "solve") {ansnam <- "inverse"} else if (minv.method == "chol2inv") {ansnam <- "inverse from Choleski decomposition"} else {ansnam <- "generalised inverse"}
    } 
    else 
    {
        ### Normal inverse
        ans <- try(solve(m), silent = TRUE)
        ansnam <- "inverse"
        if (!class(ans) == "try-error") {not.ok <- any(diag(ans) < 0) | any(eigen(ans)$values < 0)} else {not.ok <- TRUE}
        if (not.ok) 
        {
            ### Choleski inverse
            ans <- try(chol2inv(m), silent = TRUE)
            ansnam <- "inverse from Choleski decomposition"
            if (!class(ans) == "try-error") {not.ok <- any(diag(ans) < 0) | any(eigen(ans)$values < 0)} else {not.ok <- TRUE}
            if (not.ok) 
            {
                ### Generalised inverse
                ans <- try(ginv(m), silent = TRUE)
                ansnam <- "generalised inverse"
            }
        }
    }
    return(list(inv = ans, inv.nam = ansnam))
}
