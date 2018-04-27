\name{anova.jmdem}
\alias{anova.jmdem}
\title{
    Analysis of Deviance for Joint Mean and Dispersion Effect Models Fits
}
\description{
    Compute an analysis of deviance table for one or more double generalised linear model fits.
}
\usage{
\method{anova}{jmdem}(object, ..., test = NULL, type = c("1", "3"), 
      print.results = TRUE)
}
\arguments{
    \item{object, ...}{one or several objects of class \code{jmdem}, typically the result of a call to \code{jmdem}.}
    \item{test}{a character string, (partially) matching one of "\code{Rao}" or "\code{Wald}". See \code{\link{stat.anova.jmdem}}.}
    \item{type}{a character string or integer, specifying whether a type "\code{1}" (\emph{sequential}) analysis or a type "\code{3}" (\emph{partial}) analysis should be conducted. It is only relevant if a single object is specified in \code{object}. Both numeric and character inputs are allowed. See details for type 1 and type 3 analysis.}
    \item{print.results}{logical, \code{TRUE} if the result table should be printed directly, \code{FALSE} if the results should be saved in an user-defined object.}
}
\details{
    Specifying a single object gives a analysis of deviance table for that fit. If type \code{1} analysis is specified, a sequential analysis will be conducted. That is, the reductions in the residual deviance as each term of the formula is \emph{added} in turn are given in as the rows of a table, plus the residual deviances themselves.

    Type \code{3} analysis gives the reduction in the residual deviance of the fitted model after \emph{removing} each term of the formula individually, that in turn are given as the rows of a table.

    If more than one object is specified, the table has a row for the residual degrees of freedom and deviance for each model. For all but the first model, the change in degrees of freedom and deviance is also given. (This only makes statistical sense if the models are nested.) It is conventional to list the models from smallest to largest, but this is up to the user.

    The table will optionally contain "\code{Rao}" or "\code{Wald}" test statistics (and P values) comparing the model specified in the current row and the row above (type \code{1}) or the full model (type \code{3}). Both "\code{Rao}" and "\code{Wald}" test statistics are asymptotically chi-square distributed. "\code{LRT}" (Likelihood ratio test) and "\code{F}" ((F test) are not included in \code{anova.jmdem} because the comparison of the deviances of two joint mean and dispersion effects models is questionable, if not even invalid. One important argument is that the dependent variables of two different dispersion submodels given two different mean submodels are not the identical.
}
\value{
    An object of class "\code{anova}" inheriting from class "\code{data.frame}".
    
    If \code{print.results = TRUE},
    
    \item{table.x}{the \code{anova} table constructed for the mean submodel.}
    
    \item{table.z}{the \code{anova} table constructed for the dispersion submodel.}
}
\references{
    Hastie, T. J. and Pregibon, D. (1992). \emph{Generalized linear models}. Chapter 6 of Statistical Models in S eds J. M. Chambers and T. J. Hastie, Wadsworth & Brooks/Cole.

    McCullagh P. and Nelder, J.A. (1989) \emph{Generalized Linear Models}. London: Chapman and Hall.
    
    Smyth, G.K. (1989). \emph{Generalized linear models with varying dispersion}. J. R. Statist. Soc. B, \strong{51} (1), 47-60.
    
    Smyth, G.K., Verbyla, A.P. (1996). \emph{A conditional likelihood approach to residual maximum linear estimation in generalized linear models}. J. R. Statist. Soc. B, \strong{58} (3), 565-572.

    Smyth, G.K., Verbyla, A.P. (1999). \emph{Adjusted likelihood methods for modelling dispersion in generalized linear models}. Environmetrics, \strong{10}, 695-709.
    
    Wu, K.Y.K., Li, W.K. (2016). \emph{On a dispersion model with Pearson residual responses}. Comput. Statist. Data Anal., \strong{103}, 17-27.

}
\author{
    Ka Yui Karl Wu (\url{kykarlwu@gmail.com})
}
\section{Warning}{
    The comparison between two or more models will only be valid if they are fitted to the same dataset. This may be a problem if there are missing values and \code{R}'s default of  \code{na.action = na.omit} is used, and  \code{anova} will detect this with an error.
}
\seealso{
\code{\link{jmdem}}, \code{\link{anova}}
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
              
## Run a partial analysis (type 3) with Wald test
anova(fit, test = "Wald", type = 3)
}
