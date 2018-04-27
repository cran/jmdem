\name{score.jmdem, wald.jmdem}
\alias{score.jmdem}
\alias{wald.jmdem}
\title{
Asymptotic tests for fits of joint mean and dispersion effects models
}
\description{
Computes a score (Rao) or Wald chi-squared test statistics for comparing two \code{jmdem} models.
}
\usage{
score.jmdem(object, ...)

wald.jmdem(object, ...)
}
\arguments{
  \item{object}{a model or \code{list} of two or more models fitted by \code{jmdem} to be tested. Pairwise tests will be conducted.}
  
  \item{...}{a \code{list} of two or more fitted models to be tested.}
}
\details{
Given a vector of model coefficients of length p, \eqn{\Theta = (\theta_1, \ldots, \theta_q, \theta_{q+1}, \ldots, \theta_p)^T}, the score and Wald tests are usually used to test the null hypothesis against an alternative
\deqn{H_0: \theta_{q+1} = \ldots = \theta_p = 0 vs. H_0 not true}
Thus, they are asymptotic tests on the explanatory power of one or more regressors. And the result of the score and Wald tests only makes sense if the models involved are nested, i.e. all coefficients of a "smaller" (null, restricted) model are included in a "bigger" (alternative, unrestricted) model.

The main difference between the score and Wald tests is that the score test only requires the knowledge of the fitted coefficients of the "small" model. The Wald test, on the other hand, only need the estimates of the "bigger" model. Nevertheless, these tests are asymptotically equivalent, i.e. for large samples, the test statistics of these tests on the same set of models should be very close.

The key assumption is that the coefficient estimates asymptotically follow a (multivariate) normal distribution with mean and variance equal to the model parameters and their variance-covariance matrix.

\code{score.jmdem} and \code{wald.jmdem} extract the fitted coefficients and their variance-covariance matrix from the model objects, and evaluate the test statistics subsequently. So it is not necessary to specify the coefficients and variance-covariance matrix in the function arguments.

\code{score.jmdem} and \code{wald.jmdem} only return the test statistics. They are asymptotically chi-square distributed with \eqn{p - q} degrees of freedom.
}
\value{
\code{score.jmdem} and \code{wald.jmdem} return a column matrix containing the test statistics of the pairwise comparisons of all models given by the user in \code{object} and \code{...}.
}
\references{
    Engle, R.F. (1983). \emph{Wald, Likelihood Ratio, and Lagrange Multiplier Tests in Econometrics}. In Intriligator, M. D.; Griliches, Z. \emph{Handbook of Econometrics}. \strong{II}. Elsevier. pp. 796-801.

    McCullagh P. and Nelder, J.A. (1989) \emph{Generalized Linear Models}. London: Chapman and Hall.
    
    Wu, K.Y.K., Li, W.K. (2016). \emph{On a dispersion model with Pearson residual responses}. Comput. Statist. Data Anal., \strong{103}, 17-27.
}
\author{
    Ka Yui Karl Wu (\url{kykarlwu@gmail.com})
}
\note{
The score test is sometimes also called the Rao's score test or Lagrange multiplier (LM) test in different literatures.

Normally, asymptotic tests include likelihood ratio (LR), Rao's score and Wald tests. The likelihood ratio test is omitted here because the comparison of the deviances of two joint mean and dispersion effects models is questionable, if not even invalid. One important argument is that the dependent variables of two different dispersion submodels given two different mean submodels are not the identical.
}
\seealso{
\code{\link{anova.jmdem}}, \code{\link{anova}}, \code{\link{jmdem}}
}
\examples{
## Example in jmdem(...)
MyData <- simdata.jmdem.sim(mformula = y ~ x + delta, dformula = ~ z, 
                            mfamily = gaussian(), 
                            dfamily = Gamma(link = "log"), 
                            beta.true = c(0.5, 4, 1), 
                            lambda.true = c(2.5, 3), n = 100)

fit <- jmdem(mformula = y ~ x + delta, dformula = ~ z, data = MyData, 
             mfamily = gaussian, dfamily = Gamma, dev.type = "pearson", 
             method = "BFGS", df.adj = TRUE)

fit.1 <- update(fit, mformula = . ~ . - delta)

fit.2 <- update(fit.1, mformula = . ~ . - x)

## conduct a Wald tests
wald.jmdem(fit, fit.1, fit.2)

## should deliver the same results as above
wald.jmdem(object = list(fit, fit.1, fit.2))

## conduct the score test and compute the p-value directly.
raotest <- score.jmdem(fit, fit.2)
pchisq(raotest, df = abs(df.residual(fit) - df.residual(fit.2)), 
       lower.tail = FALSE)
}
