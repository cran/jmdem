\name{jmdem.summaries}
\alias{family.jmdem}
\alias{formula.jmdem}
\alias{residuals.jmdem}
\title{
    Accessing Joint Mean and Dispersion Effect Model Fits
}
\description{
    These functions are all \code{\link{methods}} for class \code{jmdem} or \code{summary.jmdem} objects.
}
\usage{
\method{formula}{jmdem}(x, submodel = c("both", "mean", "dispersion"), ...)

\method{family}{jmdem}(object, submodel = c("both", "mean", "dispersion"), ...)

\method{residuals}{jmdem}(object, type = c("deviance", "pearson", "working",
                           "response", "partial"), ...)
}
\arguments{
    \item{x, object}{the function family accesses the family objects which are stored within objects created by \code{jmdem.}}
  
    \item{submodel}{character. The family of the specified submodel. For \code{both}, the families of the mean and dispersion submodels will be return in a list of 2 elements.}
  
    \item{type}{character. For \code{residuals}, the type of residuals which should be returned. The alternatives are: "\code{deviance}" (default), "\code{pearson}", "\code{working}", "\code{response}", and "\code{partial}".}
  
    \item{...}{further arguments passed to methods.}
}
\details{
\code{family} is a generic function with methods for class "jmdem". See \code{\link{family}} for details.

Here \code{formula} is referred to the case that it is called on a fitted \code{jmdem} model object. The default first, depending on the specified \code{submodel} argument, looks for a "\code{mean.formula}" and/or "\code{dispersion.formula}" component of the \code{jmdem} object (and evaluates it), then a "\code{mean.terms}" and/or "\code{dispersion.terms}" component, then a \code{mformula} and/or \code{dformula} parameter of the call (and evaluates its value) and finally a "\code{formula}" attribute.

The references define the types of residuals: Davison & Snell is a good reference for the usages of each.

The partial residuals are a matrix of working residuals, with each column formed by omitting a term from the model.

How \code{residuals} treats cases with missing values in the original fit is determined by the \code{na.action} argument of that fit. If \code{na.action = na.omit} omitted cases will not appear in the residuals, whereas if \code{na.action = na.exclude} they will appear, with residual value \code{NA}. See also \code{\link{naresid}}.

For fits done with \code{y = FALSE} the response values are computed from other components.
}
\seealso{
\code{\link{jmdem}}, \code{\link{anova.jmdem}}, \code{\link{coef}}, \code{\link{deviance}}, \code{\link{df.residual}}, \code{\link{effects}}, \code{\link{fitted}}, \code{\link{weighted.residuals}}, \code{\link{residuals}}, \code{\link{residuals.jmdem}}, \code{\link{summary.jmdem}}, \code{\link{weights}}.
}
\references{
Cox, D. R. and Snell, E. J. (1981). \emph{Applied Statistics; Principles and Examples}. London: Chapman and Hall.

Chambers, J. M. and Hastie, T. J. (1992) \emph{Statistical Models in S}. Wadsworth & Brooks/Cole.

Davison, A. C. and Snell, E. J. (1991). \emph{Residuals and diagnostics}. In: Statistical Theory and Modelling. In Honour of Sir David Cox, FRS, eds. Hinkley, D. V., Reid, N. and Snell, E. J., Chapman & Hall.

Dobson, A. J. (1983). \emph{An Introduction to Statistical Modelling}. London: Chapman and Hall.

Hastie, T. J. and Pregibon, D. (1992). \emph{Generalized linear models}. Chapter 6 of \emph{Statistical Models in S} eds J. M. Chambers and T. J. Hastie, Wadsworth & Brooks/Cole.

McCullagh P. and Nelder, J. A. (1989). \emph{Generalized Linear Models}. London: Chapman and Hall.
}
\author{
    Karl Wu Ka Yui (karlwuky@suss.edu.sg)
}
\examples{
## The jmdem(...) example
MyData <- simdata.jmdem.sim(mformula = y ~ x, dformula = ~ z, 
                            mfamily = poisson(), 
                            dfamily = Gamma(link = "log"), 
                            beta.true = c(0.5, 4), 
                            lambda.true = c(2.5, 3), n = 100)
                            
fit <- jmdem(mformula = y ~ x, dformula = ~ z, data = MyData, 
             mfamily = poisson, dfamily = Gamma(link = "log"), 
             dev.type = "deviance", method = "CG")

coef(fit)
plot(resid(fit), fitted(fit))
abline(h = 0, lty = 2, col = "gray")
}