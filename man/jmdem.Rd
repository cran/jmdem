\name{jmdem}
\alias{jmdem}
\alias{jmdem.fit}
\title{
    Fitting Joint Mean and Dispersion Effects Models
}
\description{
    \code{jmdem} is used to fit joint mean and dispersion effects models, specified by giving a symbolic description of the linear predictors for the mean and dispersion and a description of the error distribution
}
\usage{
jmdem(mformula, dformula, data, mfamily = gaussian, dfamily = Gamma, 
      weights, subset, dev.type = c("deviance", "pearson"), 
      moffset = NULL, doffset = NULL, mustart = NULL, phistart = NULL, 
      betastart = NULL, lambdastart = NULL, hessian = TRUE, na.action, 
      grad.func = TRUE, fit.method = "jmdem.fit", 
      method = c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN", "Brent"), 
      df.adj = FALSE, disp.adj = FALSE, full.loglik = FALSE, 
      beta.first = TRUE, prefit = TRUE, mcontrasts = NULL, 
      dcontrasts = NULL, control = list(...), 
      minv.method = c("solve", "chol2inv", "ginv"), ...)

jmdem.fit(x, y, z = NULL, weights, mfamily = gaussian, dfamily = Gamma, 
          mu, phi, beta, lambda, moffset = NULL, doffset = NULL, 
          dev.type = c("deviance", "pearson"), hessian = TRUE, 
          method = c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN", "Brent"), 
          disp.adj = FALSE, df.adj = FALSE, full.loglik = FALSE, 
          control = list(), mintercept = TRUE, dintercept = TRUE, 
          grad.func = TRUE, lower = -Inf, upper = Inf, ...)
}
\arguments{
    \item{mformula}{an object of class "\code{\link{formula}}" (or one that can be coerced to that class): a symbolic description of the \emph{mean} submodel to be fitted. The details of model specification are given under 'Details'.}
  
    \item{dformula}{a symbolic description of the \emph{dispersion} submodel to be fitted. The details are also given under 'Details'.}
  
    \item{data}{an optional data frame, list or environment (or object coercible by \link{as.data.frame} to a data frame) containing the variables in the model. If not found in data, the variables are taken from \code{environment(formula)}, typically the environment from which \code{jmdem} is called.}

    \item{mfamily}{a description of the error distribution and link function to be used in the \emph{mean} submodel. This can be a character string naming a family function, a family function or the result of a call to a family function. (See \code{\link{family}} for details of family functions.)}
  
    \item{dfamily}{a description of the error distribution and link function to be used in the \emph{dispersion} submodel. (Also see \code{\link{family}} for details of family functions.)}
  
    \item{weights}{an optional vector of 'prior weights' to be used in the fitting process. Should be \code{NULL} or a numeric vector.}

    \item{subset}{an optional vector specifying a subset of observations to be used in the fitting process.}
    
    \item{dev.type}{a specification of the type of residuals to be used as the response of the \emph{dispersion} submodel. The ML estimates of the jmdem are the optima of either the quasi-likelihood function for \emph{deviance residuals}, or the pseudo-likelihood function for \emph{Pearson} residuals.}

    \item{moffset}{an a priori known component to be included in the linear predictor of the \emph{mean} submodel during fitting. This should be NULL or a numeric vector of length equal to the number of cases. One or more offset terms can be included in the formula instead or as well, and if more than one is specified their sum is used. See \code{\link{model.offset}}.}

    \item{doffset}{an a priori known component to be included in the linear predictor of the \emph{dispersion} submodel during fitting. See \code{\link{model.offset}}.}

    \item{mustart, mu}{a vector of starting values of individual means.}

    \item{phistart, phi}{a vector of starting values of individual dispersion.}

    \item{betastart, beta}{a vector of starting values for the regression parameters of the \emph{mean} submodel.}
  
    \item{lambdastart, lambda}{a vector of starting values for the regression parameters of the \emph{dispersion} submodel.}
  
    \item{hessian}{the method used to compute the information matrix. Hessian matrix will be calculated for \code{"TRUE"}, Fisher matrix for \code{"FALSE"}.}

    \item{na.action}{a function which indicates what should happen when the data contain \code{NA}s. The default is set by the \code{na.action} setting of \code{\link{options}}, and is \code{\link{na.fail}} if that is unset. The 'factory-fresh' default is \code{\link{na.omit}}. Another possible value is \code{NULL}, no action. Value \code{\link{na.exclude}} can be useful.}

    \item{grad.func}{the gradient function will be included in the optimisation for the "\code{BFGS}", "\code{CG}" and "\code{L-BFGS-B}" methods for \code{"TRUE"}. If it is \code{NULL}, a finite-difference approximation will be used.

        For the "\code{SANN}" method it specifies a function to generate a new candidate point. If it is \code{NULL} a default Gaussian Markov kernel is used.}

    \item{fit.method}{the method to be used in fitting the model. The default method "\code{jmdem.fit}" uses the general-purpose optimisation (\code{optim}): the alternative "model.frame" returns the model frame and does no fitting.

        User-supplied fitting functions can be supplied either as a function or a character string naming a function, with a function which takes the same arguments as \code{jmdem.fit}. If specified as a character string it is looked up from within the \emph{stats} namespace.}
        
    \item{method}{the method to be used for the optimisation. See \code{\link{optim}} for details.}

    \item{df.adj}{an adjustment factor for the degrees of freedom \code{(n-p)/n}, where \code{n} is the number of observations and \code{p} is the number of parameters to be estimated in \code{jmdem}, will be multiplied to the likelihood function before the optimisation for \code{"TRUE"}.}

    \item{disp.adj}{an adjustment factor for the dispersion weight will be multiplied to the estimated dispersion parameter during the optimisation for \code{"TRUE"}. For details, please see McCullagh and Nelder (1989, Ch. 10, P. 362).}

    \item{full.loglik}{the full likelihood function instead of the quasi- or pseudo-likelihood function will be used for the optimisation for \code{TRUE}.}
  
    \item{beta.first}{the mean effects will be estimated (assuming constant sample dispersion) at the initial stage for \code{TRUE}. For \code{FALSE}, the dispersion effects will be estimated first (assuming constantly zero mean for the whole sample).}

    \item{prefit}{a specfication whether \code{jmdem} uses \code{glm} to prefit the starting values of the mean and dispersion parameters. For \code{FALSE}, the initial parameter values of all the regressors are set to zero and the sample mean and sample dispersion will be used as the starting values of the corresponding submodel intercepts instead. If the submodels have no intercept, all parameters will also be set to zero. The sample mean and sample dispersion will then be used as \code{mustart} and \code{phistart} in the internal computation (they will not be officially recorded in \code{mustart} and \code{phistart} in the output object). Defaule value is \code{TRUE}.}

    \item{mcontrasts}{an optional list for the mean effect constrasts. See the \code{contrasts.arg} of \code{model.matrix.default}.}
  
    \item{dcontrasts}{an optional list for the dispersion effect constrasts. See the \code{contrasts.arg} of \code{model.matrix.default}.}
  
    \item{control}{a list of parameters for controlling the fitting process. For \code{jmdem.fit} this is passed to \code{jmdem.control}.}

    \item{minv.method}{the method used to invert matrices during the estimation process. "\code{solve}" gives the solutions of a system of equations, "\code{chol2inv}" gives the inverse from Choleski or QR decomposition and "\code{ginv}" gives the generalized inverse of a matrix. If none of the methods is specified or if they are specified in a vector such as \code{c("solve", "chol2inv", "ginv")}, the matrices will be inverted by the methods in the sequence as given in the vector until it is found.}
    
    \item{x, y, z}{\code{x} is a \emph{mean} submodel's design matrix of dimension \code{n * p}, \code{z} is a \emph{dispersion} submodel's design matrix of dimension \code{n * k}, and \code{y} is a vector of observations of length \code{n}. If \code{z} is \code{NULL}, the \emph{dispersion} submodel only contains an intercept.}

    \item{mintercept}{a specification whether the intercept term is included in the \emph{mean} submodel.}

    \item{dintercept}{a specification whether the intercept term is included in the \emph{dispersion} submodel.}
    
    \item{lower, upper}{bounds on the variables for the "\code{L-BFGS-B}" optimisation method.}

    \item{...}{For \code{control}: arguments to be used to form the default control argument if it is not supplied directly. For \code{jmdem} and \code{jmdem.fit}: further arguments passed to or from other methods.}
}
\details{
    A typical predictor has the form \code{response ~ terms} where response is the (numeric) response vector and terms is a series of terms which specifies a linear predictor for response. 
    
    A terms specification of the form \code{first + second} indicates all the terms in first together with all the terms in second with any duplicates removed. A specification of the form \code{first:second} indicates the set of terms obtained by taking the interactions of all terms in first with all terms in second. The specification \code{first * second} indicates the cross of first and second. This is the same as \code{first + second + first:second}.

    The terms in the formula will be re-ordered so that main effects come first, followed by the interactions, all second-order, all third-order and so on: to avoid this pass a terms object as the formula.

    An additional term \code{response ~ terms + eta} can be added to \code{dformula} if the \emph{mean} submodel is nested in the \emph{dispersion} submodel in the form such that \deqn{g(E(y_i))=\boldsymbol{x}_i\boldsymbol{\beta}=\eta_i, h(\phi)=\boldsymbol{z}_i\boldsymbol{\lambda}+\eta_i\gamma.} In the contrary, if the \emph{dispersion} submodel is nested in the \emph{mean} submodel such that \deqn{g(E(y_i))=\boldsymbol{x}_i\boldsymbol{\beta}+\delta_i\kappa, h(\phi_i)=\boldsymbol{z}_i\boldsymbol{\lambda}=\delta_i,} \code{mformula} can be specified as \code{response ~ terms + delta}.

    Non-\code{NULL} weights can be used to indicate that different observations have different dispersions (with the values in weights being inversely proportional to the dispersions); or equivalently, when the elements of weights are positive integers \eqn{w_i}, that each response \eqn{y_i} is the mean of \eqn{w_i} unit-weight observations. For a binomial GLM prior weights are used to give the number of trials when the response is the proportion of successes: they would rarely be used for a Poisson GLM.

    If more than one of \code{etastart} and \code{mustart} is specified, the first in the list will be used. It is often advisable to supply starting values for a quasi family, and also for families with unusual links such as gaussian("log").
    
    \code{glm.fit} is the workhorse function: it is not normally called directly but can be more efficient where the response vector, design matrix and family have already been calculated.
}
\value{
    \item{coefficients}{a named vector of estimated coefficients of both the mean and \emph{dispersion} submodel}
    \item{beta}{estimated coefficients of the \emph{mean} submodel}
    \item{lambda}{estimated coefficients of the \emph{dispersion} submodel}
    \item{residuals}{the \emph{working} residuals, that is the residuals in the final iteration of the \code{optim} fit. Depending on the type of deviance specified by \code{dev.type}, \code{residuals} corresponds to \code{deviance.residuals} or \code{pearson.residuals}. Since cases with zero weights are omitted, their working residuals are \code{NA}.}
    \item{deviance.residuals}{the \emph{deviance} residuals resulting from the final iteration of the \code{optim} fit.}
    \item{pearson.residuals}{the \emph{pearson} residuals resulting from the final iteration of the \code{optim} fit.}
    \item{fitted.values}{the fitted mean values, obtained by transforming the linear predictors by the inverse of the link function.}
    \item{dispersion}{the fitted individual dispersion values, obtained by transforming the linear predictors of the \emph{dispersion} submodel by the corresponding inverse of the link function.}
    \item{mean.rank}{the numeric rank of the fitted \emph{mean} submodel.}
    \item{dispersion.rank}{the numeric rank of the fitted \emph{dispersion} submodel.}
    \item{rank}{the total numeric rank of the fitted model. \code{mean.rank} and \code{dispersion.rank} are the corresponding ranks of the fitted \emph{mean} and \emph{dispersion} submodels.}
    \item{mean.family}{the \code{family} object used for the \emph{mean} submodel.}
    \item{dispersion.family}{the \code{family} object used for the \emph{dispersion} submodel.}
    \item{mean.linear.predictors}{the linear fit on link scale of the \emph{mean} submodel.}
    \item{dispersion.linear.predictors}{the linear fit on link scale of the \emph{dispersion} submodel.}
    \item{deviance}{the residual sum of squares of the complete fitted model.}
    \item{individual.loglik}{individual value of the log-likelihood function given the estimated mean and dispersion.}
    \item{aic}{the \emph{Akaike Information Criterion}, minus twice the maximised log-likelihood plus twice the number of parameters.}
    \item{iter}{number of iterations needed for the fit.}
    \item{weights}{the working weights, that is the weights in the final iteration of the \code{optim} fit.}
    \item{prior.weights}{the weights initially supplied, a vector of \code{1}s if none were.}
    \item{info.matrix}{the information matrix given the estimated model coefficients. The diagonal elements of its inverse are the standard errors of the model parameters.}
    \item{df.residual}{the residual degrees of freedom of the complete fitted model.}
    \item{y}{the \code{y} vector used.}
    \item{x}{the \emph{mean} submodel design matrix.}
    \item{z}{the \emph{dispersion} submodel design matrix.}
    \item{log.llh}{the maximised log-likelihood of the entire sample.}
    \item{converged}{logical. Was the \code{optim} algorithm judged to have converged?}
    \item{gradient}{logical. Was the gradient function included in the \code{optim} algorithm?}
    \item{deviance.type}{the type of redidual deviance specified, it is either "\code{deviance}" or "\code{pearson}".}
    \item{information.type}{the type of information matrix specified, it is either "\code{Hessian}" or "\code{Fisher}".}
    \item{dispersion.adjustment}{logical. Was the dispersion parameter adjusted by an adjustment factor during the optimisation?}
    \item{df.adjustment}{logical. Was the likelihood function adjusted by the degrees of freedom adjustment factor?}
    \item{optim.method}{the name of the method used in \code{optim}.}
    \item{control}{the value of the control argument used.}
    \item{data}{the evaluated dataset specified in the \code{data} argument.}
    \item{mean.model}{the model frame of the \emph{mean} submodel.}
    \item{dispersion.model}{the model frame of the \emph{dispersion} submodel.}
    \item{call}{the matched call.}
    \item{mean.formula}{the formula of the \emph{mean} submodel supplied.}
    \item{dispersion.formula}{the formula of the \emph{dispersion} submodel supplied.}
    \item{fit.method}{the name of the fit function used, currently always "\code{jmdem.fit}".}
    \item{mean.offset}{the offset vector used in the \emph{mean} submodel.}
    \item{dispersion.offset}{the offset vector used in the \emph{dispersion} submodel.}
    \item{dispersion.deviance}{the deviance sum of squares of the \emph{dispersion} submodel.}
    \item{dispersion.df.residual}{the residual degrees of freedom of the \emph{dispersion} submodel.}
    \item{null.deviance}{the residual sum of squares of the complete null model.}
    \item{df.null}{the residual degrees of freedom for the complete null model.}
    \item{dispersion.null.deviance}{the residual sum of squares of the dispersion null submodel.}
    \item{dispersion.df.null}{the residual degrees of freedom for the dispersion null submodel.}
    \item{beta.null}{the estimated coefficients of the mean null submodel.}
    \item{lambda.null}{the estimated coefficients of the dispersion null submodel.}
    \item{dispersion.null}{the estimated dispersion of the complete null model.}
    \item{residuals.null}{the residuals of the complete null model.}
    \item{mustart}{the vector of starting values for individual means used.}
    \item{phistart}{the vector of starting values for individual dispersion used.}
    \item{betastart}{the vector of starting values for the \emph{mean} submodel parameters used.}
    \item{lambdastart}{the vector of starting values for the \emph{dispersion} submodel parameters used.}
    \item{mean.terms}{the terms object used for the \emph{mean} submodel.}
    \item{dispersion.terms}{the terms object used for the \emph{dispersion} submodel.}
    \item{xlevels}{a record of the levels of the factors used in fitting the \emph{mean} submodel.}
    \item{zlevels}{a record of the levels of the factors used in fitting the \emph{dispersion} submodel.}
    \item{mean.contrasts}{the contrasts used for the \emph{mean} submodel.}
    \item{dispersion.contrasts}{the contrasts used for the \emph{dispersion} submodel.}
    \item{na.action}{information returned by model.frame on the special handling of \code{NA}s.}
    \item{init.mean.fit}{the initial values of the \emph{mean} submodel coefficients, linear predictors and fitted values.}
    \item{init.dispersion.fit}{the initial values of the \emph{dispersion} submodel coefficients, linear predictors and fitted values.}
    \item{matrix.inverse.method}{information returned on the method used for inverting matrices during optimisation.}
}

\references{
    Carroll, R.J., Ruppert, D. (1988). \emph{Transformation and Weighting in Regression}. London: Chapman and Hall.
    
    Cordeiro, M.G., Simas, A.B. (2009). \emph{The distribution of pearson residuals in generalized linear models}. Comput. Statist. Data Anal., \strong{53}, 3397-3411.
    
    %%Fletcher, R. and Reeves, C.M. (1964). \emph{Function minimization by conjugate gradients}. Computer Journal, \strong{7}, 148-154.

    McCullagh, P. (1983). \emph{Quasi-likelihood functions}. Annals of Statistics \strong{11} (1), 59-67.
    
    McCullagh P. and Nelder, J.A. (1989) \emph{Generalized Linear Models}. London: Chapman and Hall.
    
    Nash, J.C. (1990). \emph{Compact Numerical Methods for Computers. Linear Algebra and Function Minimisation}. Adam Hilger.

    Nelder, J.A., Lee, Y., Bergman, B., Hynen, A., Huele, A.F., Engel, J. (1998). \emph{Joint modelling of mean and dispersion}. Technometrics, \strong{40} (2), 168-175.
    
    %%Nelder, J.A., Mead, R. (1965). \emph{A simplex algorithm for function minimization}. Computer Journal, \strong{7}, 308-313.

    Nelder, J.A., Pregibon, D. (1987). \emph{An extended quasi-likelihood function}. Biometrika, \strong{74} (2), 221-232.
    
    Nocedal, J., Wright, S.J. (1999). \emph{Numerical Optimization}. Springer.

    Smyth, G.K. (1989). \emph{Generalized linear models with varying dispersion}. J. R. Statist. Soc. B, \strong{51} (1), 47-60.
    
    Smyth, G.K., Verbyla, A.P. (1996). \emph{A conditional likelihood approach to residual maximum linear estimation in generalized linear models}. J. R. Statist. Soc. B, \strong{58} (3), 565-572.

    Smyth, G.K., Verbyla, A.P. (1999). \emph{Adjusted likelihood methods for modelling dispersion in generalized linear models}. Environmetrics, \strong{10}, 695-709.
    
    Wedderburn, R. (1974). \emph{Quasi-likelihood functions, generalized linear models, and the Gauss-Newton method}. Biometrika, \strong{61} (3), 439-447.
    
    Wu, K.Y.K., Li, W.K. (2016). \emph{On a dispersion model with Pearson residual responses}. Comput. Statist. Data Anal., \strong{103}, 17-27.
}

\author{
    Karl Wu Ka Yui (karlwuky@suss.edu.sg)
}
\seealso{
\code{\link{anova.jmdem}}, \code{\link{summary.jmdem}}, etc. for \code{jmdem} methods, and the generic functions \code{\link{effects}}, \code{\link{fitted.values}}, and \code{\link{residuals}}.
}
\examples{
## Fit poisson counts by unnested mean and dispersion submodels. 
## Use log-links for both submodels. Set dispersion fitting based 
## on deviance residuals. Use conjugate gradient (CG) as 
## optimisation method.
MyData <- simdata.jmdem.sim(mformula = y ~ x, dformula = ~ z, 
                            mfamily = poisson(), 
                            dfamily = Gamma(link = "log"), 
                            beta.true = c(0.5, 4), 
                            lambda.true = c(2.5, 3), n = 100)
                            
fit <- jmdem(mformula = y ~ x, dformula = ~ z, data = MyData, 
             mfamily = poisson, dfamily = Gamma(link = "log"), 
             dev.type = "deviance", method = "CG")
      
## Fit Gaussian responses by nesting dispersion submodel in the mean 
## submodel. Use default link for both submodels. Set dispersion fitting 
## based on pearson residuals. Use quasi-Newton (BFGS) as optimisation 
## method. Adjust degrees of freedom for the likelihood function.
MyData <- simdata.jmdem.sim(mformula = y ~ x + delta, dformula = ~ z, 
                            mfamily = gaussian(), 
                            dfamily = Gamma(link = "log"), 
                            beta.true = c(0.5, 4, 1), 
                            lambda.true = c(2.5, 3), n = 100)

fit <- jmdem(mformula = y ~ x + delta, dformula = ~ z, data = MyData, 
             mfamily = gaussian, dfamily = Gamma, dev.type = "pearson", 
             method = "BFGS", df.adj = TRUE)
}
