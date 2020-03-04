\name{summary.jmdem.sim}
\alias{summary.jmdem.sim}
\alias{print.summary.jmdem.sim}
\title{
    Summarising JMDEM Simulations
}
\description{
    These functions are all \code{\link{methods}} for class \code{jmdem.sim} or \code{summary.jmdem.sim} objects.
}
\usage{
\method{summary}{jmdem.sim}(object, digits = max(3L, getOption("digits") - 3L), 
        scientific = FALSE, pvalue = 0.05, 
        minv.method = c("solve", "chol2inv", "ginv"), 
        other.call = FALSE, details = FALSE, ...)
              
\method{print}{summary.jmdem.sim}(x, digits = max(3L, getOption("digits") - 3L), scientific = FALSE, 
      pvalue = 0.05, signif.stars = getOption("show.signif.stars"), 
      other.call = FALSE, details = FALSE, ...)
}
\arguments{
    \item{object}{an object of class "\code{jmdem.sim}", usually, a result of a call to \code{\link{jmdem.sim}}.}
    \item{x}{an object of class "\code{summary.jmdem.sim}", usually, a result of a call to \code{summary.jmdem.sim}.}
    \item{digits}{the number of significant digits to use when printing.}
    \item{scientific}{logical; if \code{TRUE}, scientific notation is used when printing.}
    \item{pvalue}{a value between 0 and 1. It is used to compute the coverage proportion of the true parameter values by the simulated fits.}
    \item{minv.method}{the method used to invert matrices during the estimation process. "\code{solve}" gives the solutions of a system of equations, "\code{chol2inv}" gives the inverse from Choleski or QR decomposition and "\code{ginv}" gives the generalised inverse of a matrix. If none of the methods is specified or if they are specified in a vector such as \code{c("solve", "chol2inv", "ginv")}, the matrices will be inverted by the methods in the sequence as given in the vector until it is found.}
    \item{signif.stars}{logical. If \code{TRUE}, 'significance stars' are printed for each coefficient.}
    \item{other.call}{logical. If \code{true}, the rest of simulation call (i.e. without the mean and dispersion submodel formulas, families, true values) will be shown.}
    \item{details}{logical. If \code{true}, coefficients, standard errors, true parameter coverage (TRUE/FALSE) and asymptotic test statistics of each simulation will be listed.}
    \item{...}{further arguments passed to or from other methods.}
}

\details{
    The arithmetric mean of the coefficients, standard errors and coverage by the confidence intervals estimated in all simulations will be listed in a table. A detail listing of each simulation's results can be provided if required by \code{details = TRUE}. The summary also includes the averages of the Rao's score and Wald test statistics of all simulation fits.
    
    \code{print.summary.jmdem.sim} tries to be smart about formatting the coefficients, standard errors, etc according the number of significant \code{digits} (default of user-specified) or the usage of \code{scientific} notation or not.
}
\value{
    \item{digits}{the number of significant digits to use when printing.}
    \item{scientific}{logical value of using scientific notation when printing.}
    \item{details}{logical value of printing details of each simulation.}
    \item{other.call}{logical value of printing other parameters of the simulation call.}
    \item{pvalue}{numeric value between 0 and 1 used for the computation of the true parameter coverage.}
    \item{beta.true}{user-defined vector containing the true parameter values of the mean submodel.}
    \item{lambda.true}{user-defined vector containing the true parameter values of the dispersion submodel.}
    \item{simcall}{the component from \code{object}.}
    \item{mformula}{the component from \code{object}.}
    \item{dformula}{the component from \code{object}.}
    \item{mfamily}{the component from \code{object}.}
    \item{dfamily}{the component from \code{object}.}
    \item{coefficients}{mean and dispersion submodel parameter coefficients fitted in each simulation saved in a \code{data.frame}.}
    \item{stderr}{standard erros of all mean and dispersion submodel parameter coefficients estimated in each simulation saved in a \code{data.frame}.}
    \item{iterations}{a vectror containing the running numbers of each simulation.}
    \item{confint}{confidence intervals of all mean and dispersion submodel parameter coefficients estimated in each simulation saved in a \code{data.frame}.}
    \item{coverage}{the coverage of all true submodel parameters by the confidence intervals estimated in each simulation saved in a \code{data.frame}.}
    \item{asymp.test}{Rao's score and Wald test statistics of each simulation saved in a \code{data.frame}.}
    \item{average.summary}{Arithmetric means of the coefficients, standard errors, confidence interval converage estimated in all simulations saved in a \code{data.frame}.}
    \item{average.asymp.test}{(Arithmetric means of the Rao's score and Wald test statistics estimated in all simulations saved in a \code{data.frame}.}
}
\author{
    Karl Wu Ka Yui (karlwuky@suss.edu.sg)
}
\seealso{
    \code{\link{jmdem.sim}}, \code{\link{jmdem}}, \code{\link{summary}}
}
\examples{
## Example in jmdem.sim(...)
sim <- jmdem.sim(mformula = y ~ x, dformula = ~ z, beta.first = TRUE, 
                 mfamily = gaussian, dfamily = Gamma(link = "log"), 
                 x.str = list(type = "numeric", random.func = "rnorm", 
                              param = list(mean = 0, sd = 2)),
                 z.str = list(type = "factor", random.func = "runif", 
                              param = list(min = 0, max = 2)),
                 beta.true = c(1.5, 4), lambda.true = c(2.5, 3, -0.2), 
                 grad.func = TRUE, method = "BFGS", n = 50,
                 simnum = 10)
                  
## Summarise simulation
summary(sim, details = FALSE, other.call = TRUE)
}
