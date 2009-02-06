\name{dtweedie.dldphi}
\alias{dtweedie.dldphi}
\alias{dtweedie.dldphi.saddle}
\title{Tweedie Distributions}
\description{Derivatives of the log-likelihood with respect to \eqn{\phi}{phi}}
\usage{
dtweedie.dldphi(phi, mu, power, y )
dtweedie.dldphi.saddle(phi, mu, power, y )
}
\arguments{
\item{y}{vector of quantiles}
\item{mu}{the mean}
\item{phi}{the dispersion}
\item{power}{the value of \eqn{p}{power} such that the variance is 
\eqn{\mbox{var}[Y]=\phi\mu^p}{var(Y) = phi * mu^power}}
}
\value{
the value of the derivative 
\eqn{\partial\ell/\partial\phi}{d(l)/d(phi)}
where \eqn{\ell}{l} is the log-likelihood for the specified
Tweedie distribution.
\code{dtweedie.dldphi.saddle}
uses the saddlepoint approximation to determine the derivative;
\code{dtweedie.dldphi}
uses an infinite series expansion.
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
\seealso{
\code{\link{dtweedie.saddle}} 
\code{\link{dtweedie}} 
\code{\link{tweedie.profile}} 
\code{\link[statmod]{tweedie}} 
}

\examples{
### Plot dl/dphi against candidate values of phi
power <- 2
mu <- 1 
phi <- seq(2, 8, by=0.1)

set.seed(10000) # For reproducability
y <- rtweedie( 100, mu=mu, power=power, phi=3)
   # So we expect the maximum to occur at  phi=3

dldphi <- dldphi.saddle <- array( dim=length(phi))

for (i in (1:length(phi))) {
   dldphi[i] <- dtweedie.dldphi( y=y, power=power, mu=mu, phi=phi[i]) 
   dldphi.saddle[i] <- dtweedie.dldphi.saddle( y=y, power=power, mu=mu, phi=phi[i]) 
}

plot( dldphi ~ phi, lwd=2, type="l",
   ylab=expression(phi), xlab=expression(paste("dl / d",phi) ) )
lines( dldphi.saddle ~ phi, lwd=2, col=2, lty=2)
legend( "bottomright", lwd=c(2,2), lty=c(1,2), col=c(1,2),
   legend=c("'Exact' (using series)","Saddlepoint") )

# Neither are very good in this case!
   
} 
\keyword{models}

