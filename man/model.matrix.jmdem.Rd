\name{model.matrix.jmdem}
\alias{model.matrix.jmdem}
\title{
    Construct Design Matrices
}
\description{
    \code{model.matrix} creates a design (or model) matrix, e.g., by expanding factors to a set of dummy variables (depending on the contrasts) and expanding interactions similarly.
}
\usage{
\method{model.matrix}{jmdem}(object, submodel = c("both", "mean", "dispersion"), ...)
}
\arguments{
    \item{object}{the function family accesses the family objects which are stored within objects created by \code{jmdem.}}
  
    \item{submodel}{character. The family of the specified submodel. For \code{both}, the families of the mean and dispersion submodels will be return in a list of 2 elements.}
    \item{...}{further arguments passed to or from other methods.}
}
\details{
\code{model.matrix} creates a design matrix from the description given in \code{terms(object)}, using the data in \code{data} which must supply variables with the same names as would be created by a call to \code{model.frame(object)} or, more precisely, by evaluating \code{attr(terms(object), "variables")}. 
}
\value{
The design matrix for a regression-like model with the specified formula and data.

There is an attribute "\code{assign}", an integer vector with an entry for each column in the matrix giving the term in the formula which gave rise to the column. Value \code{0} corresponds to the intercept (if any), and positive values to terms in the order given by the \code{term.labels} attribute of the \code{terms} structure corresponding to object.

If there are any factors in terms in the model, there is an attribute "\code{contrasts}", a named list with an entry for each factor. This specifies the contrasts that would be used in terms in which the factor is coded by contrasts (in some terms dummy coding may be used), either as a character vector naming a function or as a numeric matrix.}

\seealso{
    \code{\link{model.frame}}, \code{\link{model.extract}}, \code{\link{terms}}

}
\references{
    Chambers, J. M. (1992). \emph{Data for models}. Chapter 3 of \emph{Statistical Models in S} eds J. M. Chambers and T. J. Hastie, Wadsworth & Brooks/Cole.
}
\author{
    Ka Yui Karl Wu (\url{kykarlwu@gmail.com})
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
              
## Extract the design matrix of the mean submodel
model.matrix(fit, submodel = "mean")
}
