\name{summary.jmdem}
\alias{summary.jmdem}
\alias{print.summary.jmdem}
\title{
    Summarising Joint Mean and Dispersion Effects Model Fits
}
\description{
    These functions are all \code{\link{methods}} for class \code{jmdem} or \code{summary.jmdem} objects.
}
\usage{
\method{summary}{jmdem}(object, correlation = FALSE, symbolic.cor = FALSE, ...)
              
\method{print}{summary.jmdem}(x, digits = max(3L, getOption("digits") - 3L), 
      scientific = FALSE, symbolic.cor = x$symbolic.cor, 
      signif.stars = getOption("show.signif.stars"), ...)
}
\arguments{
    \item{object}{an object of class "\code{jmdem}", usually, a result of a call to \code{\link{jmdem}}.}
    \item{x}{an object of class "\code{summary.jmdem}", usually, a result of a call to \code{summary.jmdem}.}
    \item{correlation}{logical; if \code{TRUE}, the correlation matrix of the estimated parameters is returned and printed.}	
    \item{digits}{the number of significant digits to use when printing.}
    \item{scientific}{logical; if \code{TRUE}, scientific notation is used when printing.}
    \item{symbolic.cor}{logical. If \code{TRUE}, print the correlations in a symbolic form (see \code{\link{symnum}}) rather than as numbers.}
    \item{signif.stars}{logical. If \code{TRUE}, 'significance stars' are printed for each coefficient.}
    \item{...}{further arguments passed to or from other methods.}
}
\details{
    \code{print.summary.jmdem} tries to be smart about formatting the coefficients, standard errors, etc. and additionally gives 'significance stars' if \code{signif.stars} is \code{TRUE}. The \code{coefficients}, \code{mean.coefficients} and \code{dispersion.coefficients} components of the result give the estimated coefficients and their estimated standard errors, together with their ratio. This third column is labelled \code{t-ratio} and a fourth column gives the two-tailed p-value corresponding to the \code{t-ratio} based on a Student t distribution.

    Aliased coefficients are omitted in the returned object but restored by the \code{print} method.

    Correlations are printed to the same decimal places specified in \code{digits} (or symbolically): to see the actual correlations print \code{summary(object)$correlation} directly.

    For more details, see \code{\link{summary.glm}}.
}
\value{
    \item{call}{the component from \code{object}.}
    \item{mean.family}{the component from \code{object}.}
    \item{dispersion.family}{the component from \code{object}.}
    \item{deviance}{the component from \code{object}.}
    \item{mean.terms}{the component from \code{object}.}
    \item{dispersion.terms}{the component from \code{object}.}
    \item{aic}{the component from \code{object}.}
    \item{mean.contrasts}{the component from \code{object}.}
    \item{dispersion.contrasts}{the component from \code{object}.}
    \item{df.residual}{the component from \code{object}.}
    \item{null.deviance}{the component from \code{object}.}
    \item{df.null}{the component from \code{object}.}
    \item{information.type}{the component from \code{object}.}
    \item{iter}{the component from \code{object}.}
    \item{mean.na.action}{the component from \code{object}.}
    \item{dispersion.na.action}{the component from \code{object}.}
    \item{deviance.resid}{the deviance residuals.}
    \item{pearson.resid}{the pearson residuals.}
    \item{resid}{the working residuals depends on the setting of \code{deviance.type}.}
    \item{coefficients}{the matrix of coefficients, standard errors, z-values and p-values. Aliased coefficients are omitted.}
    \item{mean.coefficients}{the matrix of coefficients, standard errors, z-values and p-values of the mean submodel.}
    \item{dispersion.coefficients}{the matrix of coefficients, standard errors, z-values and p-values of the dispersion submodel.}
    \item{deviance.type}{the type of redidual deviance specified, it is either "\code{deviance}" or "\code{pearson}".}
    \item{aliased}{named logical vector showing if the original coefficients are aliased.}
    \item{df}{a 3-vector of the rank of the model and the number of residual degrees of freedom, plus number of coefficients (including aliased ones).}
    \item{covariance}{the estimated covariance matrix of the estimated coefficients.}
    \item{digits}{the number of significant digits to use when printing.}
    \item{scientific}{logical value of using scientific notation when printing.}
    \item{covmat.method}{named method used to invert the covariance matrix.}
    \item{correlation}{(only if correlation is true.) The estimated correlations of the estimated coefficients.}
    \item{symbolic.cor}{(only if correlation is true.) The value of the argument symbolic.cor.}
}
\author{
    Ka Yui Karl Wu (\url{kykarlwu@gmail.com})
}
\seealso{
    \code{\link{jmdem}}, \code{\link{summary}}
}
\examples{
## Example in jmdem(...)
MyData <- simdata.jmdem.sim(mformula = y ~ x, dformula = ~ z, 
                            mfamily = poisson(), 
                            dfamily = Gamma(link = "log"), 
                            beta.true = c(0.5, 4), 
                            lambda.true = c(2.5, 3), n = 100)
                            
fit <- jmdem(mformula = y ~ x, dformula = ~ z, data = MyData, 
             mfamily = poisson, dfamily = Gamma(link = "log"), 
             dev.type = "deviance", method = "CG")
              
## Summarise fit with correlation matrix
summary(fit, correlation = TRUE, digits = 4)
}
