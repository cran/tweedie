\name{AICtweedie}
\alias{AICtweedie}
\title{Tweedie Distributions}
\description{The AIC for Tweedie models}
\usage{
	AICtweedie( glm.obj, k = 2)
}
\arguments{
	\item{glm.obj}{a fitted Tweedie \code{glm} object}
	\item{k}{numeric: the penalty per parameter to be used; the default is \eqn{k=2}}
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
	Computing the \acronym{AIC} can take a long time!
}
\author{Peter Dunn (\email{pdunn2@usc.edu.au})}
\references{
   Dunn, P. K. and Smyth, G. K. (2008).
   Evaluation of Tweedie exponential dispersion model densities by Fourier inversion.
   \emph{Statistics and Computing}, 
   \bold{18}, 73--86.
	
	Dunn, Peter K and Smyth, Gordon K (2005).
	Series evaluation of Tweedie exponential dispersion model densities
	\emph{Statistics and Computing},
	\bold{15}(4). 267--280.
	
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
test.data <- rgamma(n=200, scale=1, shape=1)

### Fit a Tweedie glm and find the AIC
m1 <- glm( test.data~1, family=tweedie(link.power=0, var.power=2) )

### A Tweedie glm with p=2 is equivalent to a gamma glm:
m2 <- glm( test.data~1, family=Gamma(link=log))

### The models are equivalent, so the AIC shoud be the same:
AICtweedie(m1)
AIC(m2)

} 
\keyword{models}

