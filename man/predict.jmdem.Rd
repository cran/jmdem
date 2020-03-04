\name{predict.jmdem}
\alias{predict.jmdem}
\title{
    Predict Method for JMDEM Fits}
\description{
    Obtains predictions and optionally estimates standard errors of those predictions from a fitted joint mean and dispersion effect model object.
}
\usage{
\method{predict}{jmdem}(object, newdata = NULL, type = c("link", "response"), 
        se.fit = FALSE, na.action = na.pass, ...)
}
\arguments{
    \item{object}{a fitted object of class inheriting from "\code{jmdem}".}

    \item{newdata}{optionally, a data frame in which to look for variables with which to predict. If omitted, the fitted linear predictors are used.}

    \item{type}{the type of prediction required. The default is on the scale of the linear predictors; the alternative "\code{response}" is on the scale of the response variable. Thus for a default binomial model the default predictions are of log-odds (probabilities on logit scale) and \code{type = "response"} gives the predicted probabilities.}

    \item{se.fit}{logical switch indicating if standard errors are required.}

    \item{na.action}{function determining what should be done with missing values in \code{newdata}. The default is to predict \code{NA}.}

    \item{...}{further arguments passed to or from other methods.}
}
\details{
    If \code{newdata} is omitted the predictions are based on the data used for the fit. In that case how cases with missing values in the original fit is determined by the \code{na.action} argument of that fit. If \code{na.action = na.omit} omitted cases will not appear in the residuals, whereas if \code{na.action = na.exclude} they will appear (in predictions and standard errors), with residual value \code{NA}. See also \code{\link{napredict}}.
}
\value{
    If \code{se.fit = FALSE}, a vector or matrix of predictions.

    If \code{se.fit = TRUE}, a list with components
    \item{fit}{Predictions, as for se.fit = FALSE.}

    \item{se.fit}{Estimated standard errors.}
}
\author{
    Karl Wu Ka Yui (karlwuky@suss.edu.sg)
}
\note{
    Variables are first looked for in \code{newdata} and then searched for in the usual way (which will include the environment of the formula used in the fit). A warning will be given if the variables found are not of the same length as those in \code{newdata} if it was supplied.
}

\seealso{
    \code{\link{jmdem}}
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
              
## Predict on the scale of the response variable with standard errors.
predict(fit, type = "response", se.fit = TRUE)

## Predict based on a new observation on the scale of the linear 
## predictors with standard errors.
predict(fit, newdata = data.frame(x = -1.5, z = 100), se.fit = TRUE)

}