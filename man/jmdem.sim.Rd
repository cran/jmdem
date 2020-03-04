\name{jmdem.sim}
\alias{jmdem.sim}
\alias{simdata.jmdem.sim}
\alias{getdata.jmdem.sim}
\title{
Simulate joint mean and dispersion effects models fits
}
\description{
    Simulate iterative \code{jmdem} fits on user-defined model settings
}
\usage{
jmdem.sim(mformula = "y ~ 1 + x", dformula = "~ 1 + z", data = NULL, 
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
          
simdata.jmdem.sim(mformula = "y ~ 1 + x", dformula = "~ 1 + z", beta.true, lambda.true, 
                   x.str = list(type = "numeric", random.func = "runif", param = list()), 
                   z.str = list(type = "numeric", random.func = "runif", param = list()), 
                   mfamily = gaussian, dfamily = Gamma, weights = NULL, n, simnum = 1, 
                   moffset = NULL, doffset = NULL)

getdata.jmdem.sim(object)
}
\arguments{
    \item{mformula}{the user-defined true mean submodel, expressed in form of an object of class "\code{\link{formula}}". The number of regressors and their interactions can be specified here, but not their true parameter values.}
  
    \item{dformula}{the user-defined true dispersion submodel. See \code{mformula}.}
  
    \item{data}{an optional data frame or list of several data frames. If no data are provided, \code{jmdem.sim} will generate its own data for simulation by \code{simdata.jmdem.sim}.}

    \item{beta.true}{a vector of the true parameter values of the mean submodel. The number of elements in \code{beta.true} must be identical with the number of parameters to be estimated in \code{mformula}, including the intercept if there exists one in the model.}
  
    \item{lambda.true}{a vector of the true parameter values of the dispersion submodel. The number of elements in \code{lambda.true} must be identical with the number of parameters to be estimated in \code{dformula}, including the intercept if there exists one in the model.}
  
    \item{mfamily}{a description of the error distribution and link function to be used in the \emph{mean} submodel. This can be a character string naming a family function, a family function or the result of a call to a family function. (See \code{\link{family}} for details of family functions.)}
  
    \item{dfamily}{a description of the error distribution and link function to be used in the \emph{dispersion} submodel. (Also see \code{\link{family}} for details of family functions.)}
  
    \item{dev.type}{a specification of the type of residuals to be used as the response of the \emph{dispersion} submodel. The ML estimates of the jmdem are the optima of either the quasi-likelihood function for \emph{deviance residuals}, or the pseudo-likelihood function for \emph{Pearson} residuals.}

    \item{x.str}{a list of user-specified structure for the generation of the mean submodel design matrix, including the \code{type} (numeric, character, logical etc.), an \code{r} function (\code{random.func}) to generate the values of the regressors and the corresponding parameters (\code{param}) to be passed on to (\code{random.func}). Note that all  parameters that belong to the same \code{random.func} must be put in a \code{list(...)}. See details.}

    \item{z.str}{a list of user-specified structure for the generation of the dispersion submodel design matrix, including the \code{type} (numeric, character, logical etc.), an \code{r} function (\code{random.func}) to generate the values of the regressors and the corresponding parameters (\code{param}) to be passed on to (\code{random.func}). Note that all  parameters that belong to the same \code{random.func} must be put in a \code{list(...)}. See details.}

    \item{n}{a numeric value specifying the sample size in each simulation.}
    
    \item{simnum}{a numeric value specifying the number of simulations.}
    
    \item{trace}{a specification whether the estimated coefficients should be printed to screen after each simulation.}

    \item{asymp.test}{a specification whether the Rao's score and Wald tests should be conducted for each simulation.}

    \item{...}{for \code{control}: arguments to be used to form the default control argument if it is not supplied directly. For \code{jmdem.sim}: further arguments passed to or from other methods.}

    The following arguments are used for JMDEM fitting. See \code{\link{jmdem}} for details.
    
    \item{weights}{an optional vector of 'prior weights' to be used in the fitting process. Should be \code{NULL} or -a numeric vector.}

    \item{moffset}{an a priori known component to be included in the linear predictor of the \emph{mean} submodel during fitting. This should be NULL or a numeric vector of length equal to the number of cases. One or more offset terms can be included in the formula instead or as well, and if more than one is specified their sum is used. See \code{\link{model.offset}}.}

    \item{doffset}{an a priori known component to be included in the linear predictor of the \emph{dispersion} submodel during fitting. See \code{\link{model.offset}}.}
    
    \item{mustart}{a vector of starting values of individual means.}

    \item{phistart}{a vector of starting values of individual dispersion.}

    \item{betastart}{a vector of starting values for the regression parameters of the \emph{mean} submodel.}
  
    \item{lambdastart}{a vector of starting values for the regression parameters of the \emph{dispersion} submodel.}
  
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
  
    \item{mcontrasts}{an optional list for the mean effect constrasts. See the \code{contrasts.arg} of \code{model.matrix.default}.}
  
    \item{dcontrasts}{an optional list for the dispersion effect constrasts. See the \code{contrasts.arg} of \code{model.matrix.default}.}
  
    \item{beta.first}{the mean effects will be estimated (assuming constant sample dispersion) at the initial stage for \code{TRUE}. For \code{FALSE}, the dispersion effects will be estimated first (assuming constantly zero mean for the whole sample).}

    \item{prefit}{a specfication whether \code{jmdem} uses \code{glm} to prefit the starting values of the mean and dispersion parameters. For \code{FALSE}, the initial parameter values of all the regressors are set to zero and the sample mean and sample dispersion will be used as the starting values of the corresponding submodel intercepts instead. If the submodels have no intercept, all parameters will also be set to zero. The sample mean and sample dispersion will then be used as \code{mustart} and \code{phistart} in the internal computation (they will not be officially recorded in \code{mustart} and \code{phistart} in the output object). Defaule value is \code{TRUE}.}

    \item{control}{a list of parameters for controlling the fitting process. For \code{jmdem.fit} this is passed to \code{jmdem.control}.}

    \item{minv.method}{the method used to invert matrices during the estimation process. "\code{solve}" gives the solutions of a system of equations, "\code{chol2inv}" gives the inverse from Choleski or QR decomposition and "\code{ginv}" gives the generalised inverse of a matrix. If none of the methods is specified or if they are specified in a vector such as \code{c("solve", "chol2inv", "ginv")}, the matrices will be inverted by the methods in the sequence as given in the vector until it is found.}
    
    \item{object}{one or several objects of class \code{jmdem.sim}, typically the result of a call to \code{jmdem.sim}.}
}
\details{
\code{jmdem.sim} simulates the fitting of datasets in which the regressors of the mean and dispersion submodels are generated according to the specification given in \code{x.str} and \code{z.str}. The response variable will be then generated according to the distribution specified in \code{mfamily} with linear predictor of the mean given by \code{mformula} and the linear predictor of the dispersion given by \code{dformula}.

The specifications in \code{x.str} and \code{z.str} are rather flexible if more than one independent variables are included in any of the submodels. For instance, if one of the two independent variables of the mean submodel is numeric generated from the normal distribution of mean 0 and standard deviation 1, and the other one is a 4-level factor \eqn{{0, 1, 2, 3}} generated from the uniform distribution, then they can be specified in a vector using \code{c(...)}, such as: \code{x.str = list(type = c("numeric", "factor"), random.func = c("rnorm", "runif"), param = c(list(mean = 0, sd = 1), list(min = 0, max = 3)))}.

Note that the higher the number of simulations specified in \code{simnum}, the more stabilised are the aggregated simulation results. The larger the sample size in each simulation, the less fluctuated are the estimated results among the simulations.

Users gain \code{simdata.jmdem.sim} higher control on the simulation by generating a number of datasets upon their own settings first, and not running \code{jmdem.sim} at the same time. By taking these steps, users also have the flexiblility to edit the datasets according their own individual requirements, before calling them in \code{jmdem.sim}.

Users can also extract the datasets used in \code{jmdem.sim} by \code{getdata.jmdem.sim}. This function is useful if the datasets are generated in \code{jmdem.sim} where users do not have access prior to the simulations. 

\code{getdata.jmdem.sim} and \code{simdata.jmdem.sim} can also be useful if the users would like to conduct various simulations with different \code{jmdem} settings on the same data.
}
\value{
An object of class \code{jmdem.sim} contains of a list of \code{jmdem} fits with full model information. That means, each element of the \code{jmdem.sim} object contains the same attributes as a \code{jmdem} object. See \emph{values} of \code{jmdem} for details.
}
\author{
    Karl Wu Ka Yui (karlwuky@suss.edu.sg)
}
\seealso{
    \code{\link{jmdem}}, \code{\link{summary.jmdem.sim}}
}
\examples{
## Run 10 JMDEM simulations with samples of size 50. The response
## variable is Gaussian with mean beta_0 + beta_1 * x and variance 
## log(sigma^2) = lambda_0 + lambda_1 * z. The observations of 
## the predictor x should be random numbers generated from the normal 
## distribution with mean 0 and standard deviation 2. The observations
## of z are factors with three levels between 0 and 2, generated from 
## the uniform distribution. The true values of the mean submodel's 
## intercept and slope are 1.5 and 4, as well as 2.5, 3 and -0.2 for 
## the dispersion submodel's intercept and slope.
sim <- jmdem.sim(mformula = y ~ x, dformula = ~ z, beta.first = TRUE, 
                 mfamily = gaussian, dfamily = Gamma(link = "log"), 
                 x.str = list(type = "numeric", random.func = "rnorm", 
                              param = list(mean = 0, sd = 2)),
                 z.str = list(type = "factor", random.func = "runif", 
                              param = list(min = 0, max = 2)),
                 beta.true = c(1.5, 4), lambda.true = c(2.5, 3, -0.2), 
                 grad.func = TRUE, method = "BFGS", n = 50,
                 simnum = 10)
}
