\name{tweedie.convert}
\alias{tweedie.convert}
\title{Convert Tweedie parameters}
\description{Converts Tweedie distribution parameters to the parameters of the underlying distributions}
\usage{
	tweedie.convert( xi=NULL, mu, phi, power=NULL)
}
\arguments{
	\item{xi}{the value of \eqn{\xi}{xi} such that the variance is 
	\eqn{\mbox{var}[Y]=\phi\mu^{\xi}}{var(Y) = phi * mu^xi}}
   \item{power}{a synonym for \eqn{\xi}{xi}}
	\item{mu}{the mean}
	\item{phi}{the dispersion}
}
\value{
  a list containing the values of
  the mean of the underlying Poisson distribution (as \code{poisson.lambda}),
  the scale parameter of the underlying gamma distribution (as \code{gamma.scale}),
  the shape parameter of the underlying gamma distribution (as \code{gamma.shape}),
  the probability of obtaining a zero response (as \code{p0}),
  the mean of the underlying gamma distribution (as \code{gamma.mean}),
  and
  the dispersion parameter of the underlying gamma distribution (as \code{gamma.phi}).
}
\details{
    The Tweedie family of distributions with \eqn{1<\xi<2}{1 < power < 2}
    is the Poisson sum of gamma distributions
    (where the Poisson distribution has mean \eqn{\lambda}{lambda},
    and the gamma distribution has scale and shape parameters).
    When used to fit a glm,    
    the model is fitted with the usual glm parameters:
    the mean \eqn{\mu}{mu}
    and the dispersion parameter \eqn{\phi}{phi}.
    This function converts the parameters
    \eqn{(p, \mu, \phi)}{(p, mu, phi)}
    to the values of the parameters of the underlying Poisson distribution \eqn{\lambda}{lambda}
    and gamma distribution (scale and shape parameters).
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
	
  Tweedie, M. C. K. (1984).
  An index which distinguishes between some important exponential families.
  \emph{Statistics: Applications and New Directions.
  Proceedings of the Indian Statistical Institute Golden Jubilee International Conference}
  (Eds. J. K. Ghosh and J. Roy), pp. 579-604. Calcutta: Indian Statistical Institute.
}
\seealso{\code{\link{dtweedie.saddle}} }

\examples{
tweedie.convert(xi = 1.5, mu = 1, phi = 1)
} 
\keyword{models}

