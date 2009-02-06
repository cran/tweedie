\name{tweedie.plot}
\alias{tweedie.plot}
\title{Tweedie Distributions: plotting}
\description{Plotting Tweedie density and distribution functions}
\usage{
   tweedie.plot(y, power, mu, phi, type="pdf", add=FALSE, ...)
}
\arguments{
   \item{y}{vector of values at which to evaluate and plot}
   \item{power}{the value of \eqn{p}{power} such that the variance is 
      \eqn{\mbox{var}[Y]=\phi\mu^p}{var(Y) = phi * mu^power}}
   \item{mu}{the mean}
   \item{phi}{the dispersion}
   \item{type}{what to plot: \code{pdf} (the default) means the probability function,
               or \code{cdf}, the cumulative distribution function}
   \item{add}{if \code{TRUE}, the plot is added to the current device;
              if \code{FALSE} (the default), a new plot is produced}
   \item{...}{Arguments to be passed to the plotting method}
}
\value{
   this function is usually called for side-effect of
   producing a plot of the specified Tweedie distribution,
   properly plotting the exact zero that occurs at \eqn{y=0}{y=0}
   when \eqn{1<p<2}{1<p<2}.
   However,
   it also produces a list with the computed density at the given points,
   with components \code{y} and \code{x} respectively,
   such that \code{plot(y~x)} approximately reproduces the plot.
}
\details{
   For details, see \code{\link{dtweedie}}
}
\author{Peter Dunn (\email{pdunn2@usc.edu.au})}
\references{
   Dunn, P. K. \& Smyth, G. K. (2008).
   Evaluation of Tweedie exponential dispersion model densities by Fourier inversion.
   \emph{Statistics and Computing}, 
   \bold{18}, 73--86.
	
   Dunn, Peter K and Smyth, Gordon K (2005).
   Series evaluation of Tweedie exponential dispersion model densities
   \emph{Statistics and Computing},
   \bold{15}(4). 267--280.
   
   Dunn, Peter K and Smyth, Gordon K (2001).
   Tweedie family densities: methods of evaluation.
   \emph{Proceedings of the 16th International Workshop on Statistical Modelling},
   Odense, Denmark, 2--6 July
   
   Jorgensen, B. (1987).
   Exponential dispersion models.
   \emph{Journal of the Royal Statistical Society}, B,
   \bold{49}, 127--162.
   
   Jorgensen, B. (1997).
   \emph{Theory of Dispersion Models}.
   Chapman and Hall, London.
   
   Nolan, John P (1997).
   Numerical calculation of stable densities and distribution functions.
   \emph{Communication in Statistics---Stochastic models},
   \bold{13}(4). 759--774.
   
   Sidi, Avram (1982).
   The numerical evaluation of very oscillatory infinite integrals by
   extrapolation.
   \emph{Mathematics of Computation}
   \bold{38}(158), 517--529.
   
   Sidi, Avram (1988).  
   A user-friendly extrapolation method for
   oscillatory infinite integrals.
   \emph{Mathematics of Computation}
   \bold{51}(183), 249--266.
   
   Tweedie, M. C. K. (1984).
   An index which distinguishes between some important exponential families.
   \emph{Statistics: Applications and New Directions.
   Proceedings of the Indian Statistical Institute Golden Jubilee International Conference}
   (Eds. J. K. Ghosh and J. Roy), pp. 579-604. Calcutta: Indian Statistical Institute.
}
\seealso{\code{\link{dtweedie}} }

\examples{
### Plot a Tweedie density with 1<p<2
yy <- seq(0,5,length=100)
tweedie.plot( power=1.7, mu=1, phi=1, y=yy, lwd=2)
tweedie.plot( power=1.2, mu=1, phi=1, y=yy, add=TRUE, lwd=2, col="red")
legend("topright",lwd=c(2,2), col=c("black","red"), pch=c(19,19),
   legend=c("p=1.7","p=1.2") )

### Plot distribution functions
tweedie.plot( power=1.05, mu=1, phi=1, y=yy, lwd=2, type="cdf", ylim=c(0,1))
tweedie.plot( power=2, mu=1, phi=1, y=yy, add=TRUE, lwd=2, type="cdf",col="red")
legend("bottomright",lwd=c(2,2), col=c("black","red"),
   legend=c("p=1.05","p=2") )

### Now, plot two densities, combining p>2 and 1<p<2
tweedie.plot( power=3.5, mu=1, phi=1, y=yy, lwd=2)
tweedie.plot( power=1.5, mu=1, phi=1, y=yy, lwd=2, col="red", add=TRUE)
legend("topright",lwd=c(2,2), col=c("black","red"), pch=c(NA,19),
   legend=c("p=3.5","p=1.5") )
} 
\keyword{models}
