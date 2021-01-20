\name{logLiktweedie}
\alias{logLiktweedie}
\title{Tweedie Distributions}
\description{The log likelihood for Tweedie models}
\usage{
	logLiktweedie( glm.obj, dispersion=NULL)
}
\arguments{
	\item{glm.obj}{a fitted Tweedie \code{glm} object}
	\item{dispersion}{the dispersion parameter \eqn{\phi}{phi}; the default is \code{NULL} which means to use an estimate}
}
\value{
	Returns the log-likelihood from the specified model
}
\details{
	The log-likelihood is computed from the \acronym{AIC},
	so see \code{\link{AICtweedie}} for more details.
}
    
\section{Note}{
	Computing the log-likelihood may take a long time.
}
\section{Note}{
	Tweedie distributions with the index parameter as 1
	correspond to Poisson distributions when \eqn{\phi = 1}{phi=1}.
	However,
	in general a Tweedie distribution with an index parameter equal to one
	may not be referring to a Poisson distribution with \eqn{\phi=1}{phi=1},
	so we cannot assume that \eqn{\phi=1}{phi=1} just because the index parameter is set to one.
	If the Poisson distribution is intended,
	then \code{dispersion=1} should be specified.
	The same argument applies for similar situations.
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

	Jorgensen, B. (1997).
	\emph{Theory of Dispersion Models}.
	Chapman and Hall, London.
	
   Sakamoto, Y., Ishiguro, M., and Kitagawa G. (1986). 
   \emph{Akaike Information Criterion Statistics}. 
   D. Reidel Publishing Company.
}
\seealso{\code{\link{AICtweedie}} }

\examples{
library(statmod) # Needed to use  tweedie  family object

### Generate some fictitious data
test.data <- rgamma(n = 200, scale = 1, shape = 1)

### Fit a Tweedie glm and find the AIC
m1 <- glm( test.data ~ 1, family = tweedie(link.power = 0, var.power = 2) )

### A Tweedie glm with p=2 is equivalent to a gamma glm:
m2 <- glm( test.data ~ 1, family = Gamma(link = log))

### The models are equivalent, so the AIC shoud be the same:
logLiktweedie(m1)
logLik(m2)

} 
\keyword{models}

