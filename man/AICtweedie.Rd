\name{AICtweedie}
\alias{AICtweedie}
\title{Tweedie Distributions}
\description{The AIC for Tweedie glms}
\usage{
	AICtweedie( glm.obj, dispersion=NULL, k = 2, verbose=TRUE)
}
\arguments{
	\item{glm.obj}{a fitted Tweedie \code{glm} object}
	\item{dispersion}{the dispersion parameter \eqn{\phi}{phi}; the default is \code{NULL} which means to use an estimate}
	\item{k}{numeric: the penalty per parameter to be used; the default is \eqn{k=2}}
	\item{verbose}{if \code{TRUE} (the default), a warning message is produced about the Poisson case; see the second Note below}
}
\value{
	Returns a numeric value with the
   corresponding \acronym{AIC} (or \acronym{BIC}, depending on \eqn{k})
}
\details{
	See \code{\link{AIC}} for more details on the \acronym{AIC};
	see \code{\link{dtweedie}} for more details on computing the Tweedie densities
}
    
\section{Note}{
	Computing the \acronym{AIC} may take a long time.
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
\seealso{\code{\link{AIC}} }

\examples{
library(statmod) # Needed to use  tweedie  family object

### Generate some fictitious data
test.data <- rgamma(n = 200, scale = 1, shape = 1)

### Fit a Tweedie glm and find the AIC
m1 <- glm( test.data ~ 1, family = tweedie(link.power = 0, var.power = 2) )

### A Tweedie glm with p=2 is equivalent to a gamma glm:
m2 <- glm( test.data ~ 1, family = Gamma(link = log))

### The models are equivalent, so the AIC shoud be the same:
AICtweedie(m1)
AIC(m2)

} 
\keyword{models}

