\name{stat.anova.jmdem}
\alias{stat.anova.jmdem}
\title{
    JMDEM Anova Statistics
}
\description{
    This is a utility function, used in \code{jmdem} method for \code{\link{anova}(..., test != NULL)} and should not be used by the average user.
}
\usage{
stat.anova.jmdem(table, test = c("Rao", "Wald")) 
}
\arguments{
  \item{table}{numeric matrix as results from \code{\link{anova.jmdem}(..., test = NULL, print.results = FALSE)} saved as the attributes \code{table.x} or \code{table.z}.}
  
  \item{test}{a character string, partially matching one of "\code{Rao}" or "\code{Wald}".}
}
\value{
    A matrix which is the original \code{table}, augmented by a column of test statistics, depending on the \code{test} argument.
}
\references{
    Hastie, T. J. and Pregibon, D. (1992). \emph{Generalized linear models}. Chapter 6 of \emph{Statistical Models in S} eds J. M. Chambers and T. J. Hastie, Wadsworth & Brooks/Cole.
}
\author{
    Ka Yui Karl Wu (\url{kykarlwu@gmail.com})
}
\seealso{
    \code{\link{anova.jmdem}}
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
              
stat.anova.jmdem(anova(fit, test = "Rao", print.results = FALSE)$table.x)
}
