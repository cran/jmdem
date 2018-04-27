\name{jmdem.control}
\alias{jmdem.control}
\title{
    Auxiliary for Controlling JMDEM Fitting
}
\description{
    Auxiliary function for \code{\link{jmdem}} fitting. Typically only used internally by \code{\link{jmdem.fit}}, but may be used to construct a \code{control} argument to either function.
}
\usage{
jmdem.control(maxit = 100, epsilon = 1e-8, prefit.trace = FALSE, 
              fit.trace = FALSE, null.approx = 1e-8, trace = 0, 
              fnscale = -1, parscale = 1, ndeps = 0.001, 
              abstol = -Inf, reltol = sqrt(.Machine$double.eps), 
              alpha = 1, beta = 0.5, gamma = 2, REPORT = 10, 
              type = 1, lmm = 5, factr = 1e+07, pgtol = 0, 
              temp = 10, tmax = 10)
}
\arguments{
    \item{maxit}{integer giving the maximal number of optimisation iterations.}
    \item{epsilon}{positive convergence tolerance \eqn{\epsilon}; the iterations converge when \eqn{|dev - dev_{old}|/(|dev| + 0.1) < \epsilon}.}
    \item{prefit.trace}{logical indicating if output should be produced for each iteration in the \code{prefit} process.}
    \item{fit.trace}{logical indicating if output should be produced for each iteration in the \code{jmdem.fit} process.}
    \item{null.approx}{approximisation of zeros to avoid estimation abortion in the case of \eqn{log(0)} or \eqn{1/0}.}
    The following control arguments are used by \code{optim}. Please refer to \code{\link{optim}} for details
    \item{trace}{non-negative integer. If positive, tracing information on the progress of the optimisation is produced. Higher values may produce more tracing information: for method "\code{L-BFGS-B}" there are six levels of tracing.}
    \item{fnscale}{An overall scaling to be applied to the value of \code{fn} and \code{gr} during optimisation. If negative, turns the problem into a maximisation problem. Optimisation is performed on \code{fn(par)/fnscale}.}
    \item{parscale}{A vector of scaling values for the parameters. Optimisation is performed on \code{par/parscale} and these should be comparable in the sense that a unit change in any element produces about a unit change in the scaled value. Not used (nor needed) for \code{method = "Brent"}.}
    \item{ndeps}{A vector of step sizes for the finite-difference approximation to the gradient, on \code{par/parscale} scale. Defaults to \code{1e-3}.}
    \item{abstol}{The absolute convergence tolerance. Only useful for non-negative functions, as a tolerance for reaching zero.}
    \item{reltol}{Relative convergence tolerance. The algorithm stops if it is unable to reduce the value by a factor of \code{reltol * (abs(val) + reltol)} at a step. Defaults to \code{sqrt(.Machine$double.eps)}, typically about \code{1e-8}.} 
    \item{alpha, beta, gamma}{Scaling parameters for the "\code{Nelder-Mead}" method. \code{alpha} is the reflection factor (default 1.0), \code{beta} the contraction factor (0.5) and \code{gamma} the expansion factor (2.0).}
    \item{REPORT}{The frequency of reports for the "\code{BFGS}", "\code{L-BFGS-B}" and "\code{SANN}" methods if \code{control$trace} is positive. Defaults to every 10 iterations for "\code{BFGS}" and "\code{L-BFGS-B}", or every 100 temperatures for "\code{SANN}".}
    \item{type}{for the conjugate-gradients ("\code{CG}") method. Takes value \code{1} for the Fletcher-Reeves update, \code{2} for Polak-Ribiere and \code{3} for Beale-Sorenson.}
    \item{lmm}{is an integer giving the number of \code{BFGS} updates retained in the "\code{L-BFGS-B}" method, It defaults to \code{5}.}
    \item{factr}{controls the convergence of the "\code{L-BFGS-B}" method. Convergence occurs when the reduction in the objective is within this factor of the machine tolerance. Default is \code{1e7}, that is a tolerance of about \code{1e-8}.}
    \item{pgtol}{helps control the convergence of the "\code{L-BFGS-B}" method. It is a tolerance on the projected gradient in the current search direction. This defaults to zero, when the check is suppressed.} 
    \item{tmax}{controls the "\code{SANN}" method. It is the starting temperature for the cooling schedule. Defaults to \code{10}.}
    \item{temp}{is the number of function evaluations at each temperature for the "\code{SANN}" method. Defaults to \code{10}.}
}

\details{
    The control argument of \code{jmdem} is by default passed to the control argument of {jmdem.fit}, which uses its elements as arguments to {jmdem.control}: the latter provides defaults and sanity checking.

    When \code{trace} is true, calls to \code{\link{cat}} produce the output for each iteration. Hence, \code{options(digits = *)} can be used to increase the precision, see the example.}

\value{
    A list with components named as the arguments.
}
\references{
    Belisle, C.J.P. (1992). \emph{Convergence theorems for a class of simulated annealing algorithms on Rd}. Journal of Applied Probability, \strong{29}, 885-895.

    Byrd, R. H., Lu, P., Nocedal, J. and Zhu, C. (1995). \emph{A limited memory algorithm for bound constrained optimisation}. SIAM Journal on Scientific Computing, \strong{16}, 1190-1208.

    Fletcher, R. and Reeves, C.M. (1964). \emph{Function minimization by conjugate gradients}. Computer Journal, \strong{7}, 148-154.

    Hastie, T. J. and Pregibon, D. (1992). \emph{Generalized linear models}. Chapter 6 of Statistical Models in S eds J. M. Chambers and T. J. Hastie, Wadsworth & Brooks/Cole.

    Nash, J.C. (1990). \emph{Compact Numerical Methods for Computers. Linear Algebra and Function Minimisation}. Adam Hilger.

    Nelder, J.A., Mead, R. (1965). \emph{A simplex algorithm for function minimization}. Computer Journal, \strong{7}, 308-313.

    Nocedal, J., Wright, S.J. (1999). \emph{Numerical Optimisation}. Springer.

    Smyth, G.K. (1989). \emph{Generalised linear models with varying dispersion}. J. R. Statist. Soc. B, \strong{51} (1), 47-60.
    
    Smyth, G.K., Verbyla, A.P. (1999). \emph{Adjusted likelihood methods for modelling dispersion in generalised linear models}. Environmetrics, \strong{10}, 695-709.
    
    Wu, K.Y.K., Li, W.K. (2016). \emph{On a dispersion model with Pearson residual responses}. Comput. Statist. Data Anal., \strong{103}, 17-27.
}
\author{
    Ka Yui Karl Wu (\url{kykarlwu@gmail.com})
}
\seealso{
    \code{\link{jmdem.fit}}, the fitting procedure used by \code{\link{jmdem}}.
}
\examples{
## Example in jmdem(...). Limit maximum iteration number to 20 and 
## trace the deviance development in the fitting process
MyData <- simdata.jmdem.sim(mformula = y ~ x, dformula = ~ s, 
                            mfamily = poisson(), 
                            dfamily = Gamma(link = "log"), 
                            beta.true = c(0.5, 4), 
                            lambda.true = c(2.5, 3), n = 100)
                            
fit <- jmdem(mformula = y ~ x, dformula = ~ s, data = MyData, 
             mfamily = poisson, dfamily = Gamma(link = "log"), 
             dev.type = "deviance", method = "CG",
             control = list(maxit = 20, fit.trace = TRUE))

## Change to a small convergence tolerance and trace the optimisation 
## process in optim
jmdem.control(list(epsilon = 1e-14, trace = 1))
}