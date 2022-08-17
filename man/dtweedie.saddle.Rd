\name{dtweedie.saddle}
\alias{dtweedie.saddle}
\title{Tweedie Distributions (saddlepoint approximation)}
\description{Saddlepoint density for the Tweedie distributions}
\usage{dtweedie.saddle(y, xi=NULL, mu, phi, eps=1/6, power=NULL)}
\arguments{
	\item{y}{the vector of responses}
	\item{xi}{the value of \eqn{\xi}{xi} such that the variance is 
	\eqn{\mbox{var}[Y]=\phi\mu^{\xi}}{var(Y) = phi * mu^xi}}
   \item{power}{a synonym for \eqn{\xi}{xi}}
	\item{mu}{the mean}
	\item{phi}{the dispersion}
	\item{eps}{the offset in computing the variance function.
	The default is \code{eps=1/6}
	(as suggested by Nelder and Pregibon, 1987).}
	}
	\value{saddlepoint (approximate) density
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
    
    When \eqn{1<p<2}{1 < power < 2},
    the distribution are continuous for \eqn{Y}{Y} greater than zero,
    with a positive mass at \eqn{Y=0}{Y=0}.
    For \eqn{p>2}{power > 2}, 
    the distributions are continuous for \eqn{Y}{Y} greater than zero.

    This function approximates the density using the
    saddlepoint approximation defined by Nelder and Pregibon (1987).
}
\author{Peter Dunn (\email{pdunn2@usc.edu.au})}
\references{
	Daniels, H. E. (1954).
	Saddlepoint approximations in statistics.
	\emph{Annals of Mathematical Statistics},
	\bold{25}(4), 631--650.
	
	Daniels, H. E. (1980).
	Exact saddlepoint approximations.
	\emph{Biometrika},
	\bold{67}, 59--63.
	\doi{10.1093/biomet/67.1.59}
	
   Dunn, P. K. and Smyth, G. K. (2008).
   Evaluation of Tweedie exponential dispersion model densities by Fourier inversion.
   \emph{Statistics and Computing}, 
   \bold{18}, 73--86.
	\doi{10.1007/s11222-007-9039-6}
	
	Dunn, Peter K and Smyth, Gordon K (2001).
	Tweedie family densities: methods of evaluation.
	\emph{Proceedings of the 16th International Workshop on Statistical Modelling},
	Odense, Denmark, 2--6 July
	
	Dunn, Peter K and Smyth, Gordon K (2005).
	Series evaluation of Tweedie exponential dispersion model densities
	\emph{Statistics and Computing},
	\bold{15}(4). 267--280.
	\doi{10.1007/s11222-005-4070-y}

	Jorgensen, B. (1987).
	Exponential dispersion models.
	\emph{Journal of the Royal Statistical Society}, B,
	\bold{49}, 127-162.
	
	Jorgensen, B. (1997).
	\emph{Theory of Dispersion Models},
	Chapman and Hall, London.
	
	Nelder, J. A. and Pregibon, D. (1987).
	An extended quasi-likelihood function.
	\emph{Biometrika},
	\bold{74}(2), 221--232.
	\doi{10.1093/biomet/74.2.221}
	
	Tweedie, M. C. K. (1984). 
	An index which distinguishes between some important exponential families.
	\emph{Statistics: Applications and New Directions.
	Proceedings of the Indian Statistical Institute Golden Jubilee International Conference}
	(Eds. J. K. Ghosh and J. Roy), pp. 579-604. Calcutta: Indian Statistical Institute.
}
\seealso{\code{\link{dtweedie}} }

\examples{
p <- 2.5
mu <- 1
phi <- 1
y <- seq(0, 10, length=100)
fy <- dtweedie( y=y, power=p, mu=mu, phi=phi)
plot(y, fy, type="l")
# Compare to the saddlepoint density
f.saddle <- dtweedie.saddle( y=y, power=p, mu=mu, phi=phi)
lines( y, f.saddle, col=2 )
}
\keyword{models}
