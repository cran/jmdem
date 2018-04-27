\name{update.jmdem}
\alias{update.jmdem}
\title{
    Update and Re-fit a JMDEM Call
}
\description{
    \code{update} will update and (by default) re-fit a model. It does this by extracting the call stored in the object, updating the call and (by default) evaluating that call. Sometimes it is useful to call update with only one argument, for example if the data frame has been corrected.
}
\usage{
\method{update}{jmdem}(object, mformula, dformula, ...)
}
\arguments{
    \item{object}{An existing fit from a \code{jmdem} model function}
    
    \item{mformula}{Changes to the formula of the mean submodel - see update.formula for details.}

    \item{dformula}{Changes to the formula of the dispersion submodel - see update.formula for details.}
    
    \item{...}{Additional arguments to the call, or arguments with changed values. Use \code{name = NULL} to remove the argument name.}
}
\references{
    Chambers, J. M. (1992). \emph{Linear models}. Chapter 4 of \emph{Statistical Models} in S eds J. M. Chambers and T. J. Hastie, Wadsworth & Brooks/Cole.
}
\author{
    Ka Yui Karl Wu (\url{kykarlwu@gmail.com})
}
\seealso{
    \code{\link{update.formula}}
}
\examples{
## Example in jmdem(...): Update the dispersion fitting based on Pearson
## residuals and change from Nelder-Mead to BFGS as optimisation method.
MyData <- simdata.jmdem.sim(mformula = y ~ x, dformula = ~ z, 
                            mfamily = poisson(), 
                            dfamily = Gamma(link = "log"), 
                            beta.true = c(0.5, 4), 
                            lambda.true = c(2.5, 3), n = 100)
                            
fit <- jmdem(mformula = y ~ x, dformula = ~ z, data = MyData, 
             mfamily = poisson, dfamily = Gamma(link = "log"), 
             dev.type = "deviance", method = "CG")
             
update(fit, dev.type = "pearson", method = "BFGS")
}
