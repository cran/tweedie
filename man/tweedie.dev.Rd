\name{tweedie.dev}
\alias{tweedie.dev}
\title{Tweedie Distributions: the deviance function}
\description{The deviance function for the Tweedie family of distributions}
\usage{
tweedie.dev(y, mu, power)
}
\arguments{
\item{y}{vector of quantiles (which can be zero if \eqn{1<p<2}{1 < power < 2}}
\item{mu}{the mean}
\item{power}{the value of \eqn{p}{power} such that the variance is 
\eqn{\mbox{var}[Y]=\phi\mu^p}{var(Y) = phi * mu^power}}
}
\value{
the value of the deviance 
for the given Tweedie distribution with parameters 
\code{mu}, 
\code{phi}
and 
\code{power}.
}
\details{
    The Tweedie family of distributions belong to the class
    of exponential dispersion models (\acronym{EDM}s),
    famous for their role in generalized linear models.
    The Tweedie distributions are the \acronym{EDM}s with a variance of the form
    \eqn{\mbox{var}[Y]=\phi\mu^p}{var(Y) = phi * mu^power}
    where \eqn{p}{power} is greater than or equal to one, or less than or equal to zero.
    \bold{This function only evaluates for \eqn{p}{power}
    greater than or equal to one.}
    Special cases include the 
    normal (\eqn{p=0}{power=0}),
    Poisson (\eqn{p=1}{power=1} with \eqn{\phi=1}{phi=1}),
    gamma (\eqn{p=2}{power=2})
    and
    inverse Gaussian (\eqn{p=3}{power=3})
    distributions.
    For other values of  \code{power},
    the distributions are still defined but cannot be written in closed form,
    and hence evaluation is very difficult.
    
    The deviance is defined by \code{\link{deviance}} as
    ``up to a constant, minus twice the maximized log-likelihood.
    Where sensible, the constant is chosen so that a saturated
    model has deviance zero.''

}
\author{Peter Dunn (\email{pdunn2@usc.edu.au})}
\references{
   Dunn, P. K. and Smyth, G. K. (2008).
   Evaluation of Tweedie exponential dispersion model densities by Fourier inversion.
   \emph{Statistics and Computing}, 
   \bold{18}, 73--86.
   \doi{10.1007/s11222-007-9039-6}
	
	Dunn, Peter K and Smyth, Gordon K (2005).
	Series evaluation of Tweedie exponential dispersion model densities
	\emph{Statistics and Computing},
	\bold{15}(4). 267--280.
	\doi{10.1007/s11222-005-4070-y}
	
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
	
	Sidi, Avram (1982).
	The numerical evaluation of very oscillatory infinite integrals by
	extrapolation.
	\emph{Mathematics of Computation}
	\bold{38}(158), 517--529.
	\doi{10.1090/S0025-5718-1982-0645667-5}
	
	Sidi, Avram (1988).  
	A user-friendly extrapolation method for
	oscillatory infinite integrals.
	\emph{Mathematics of Computation}
	\bold{51}(183), 249--266.
	\doi{10.1090/S0025-5718-1988-0942153-5}
	
	Tweedie, M. C. K. (1984).
	An index which distinguishes between some important exponential families.
	\emph{Statistics: Applications and New Directions.
	Proceedings of the Indian Statistical Institute Golden Jubilee International Conference}
	(Eds. J. K. Ghosh and J. Roy), pp. 579-604. Calcutta: Indian Statistical Institute.
}
\seealso{
\code{\link{dtweedie}},
\code{\link{dtweedie.saddle}},
\code{\link[statmod]{tweedie}},
\code{\link{deviance}},
\code{\link{glm}} 
}

\examples{
### Plot a Tweedie deviance function when 1<p<2
mu <- 1 

y <- seq(0, 6, length=100) 

dev1 <- tweedie.dev( y=y, mu=mu, power=1.1) 
dev2 <- tweedie.dev( y=y, mu=mu, power=1.5)
dev3 <- tweedie.dev( y=y, mu=mu, power=1.9) 

plot(range(y), range( c(dev1, dev2, dev3)), 
   type="n", lwd=2, ylab="Deviance", xlab=expression(italic(y)) )

lines( y, dev1, lty=1, col=1, lwd=2 )
lines( y, dev2, lty=2, col=2, lwd=2 )
lines( y, dev3, lty=3, col=3, lwd=2 )


legend("top", col=c(1,2,3), lwd=c(2,2,2), lty=c(1,2,3),
    legend=c("p=1.1","p=1.5", "p=1.9") )


### Plot a Tweedie deviance function when p>2
mu <- 1 

y <- seq(0.1, 6, length=100) 

dev1 <- tweedie.dev( y=y, mu=mu, power=2) # Gamma
dev2 <- tweedie.dev( y=y, mu=mu, power=3) # Inverse Gaussian
dev3 <- tweedie.dev( y=y, mu=mu, power=4) 

plot(range(y), range( c(dev1, dev2, dev3)), 
   type="n", lwd=2, ylab="Deviance", xlab=expression(italic(y)) )

lines( y, dev1, lty=1, col=1, lwd=2 )
lines( y, dev2, lty=2, col=2, lwd=2 )
lines( y, dev3, lty=3, col=3, lwd=2 )


legend("top", col=c(1,2,3), lwd=c(2,2,2), lty=c(1,2,3),
    legend=c("p=2 (gamma)", "p=3 (inverse Gaussian)", "p=4") )

} 

\keyword{models}

