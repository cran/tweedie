\name{tweedie-package}
\alias{tweedie-package}
\alias{Tweedie}
\docType{package}
\title{
Tweedie Distributions
}
\description{
Functions for computing and fitting the Tweedie family of distributions
}
\details{
\tabular{ll}{
	Package: \tab tweedie\cr
	Type: \tab Package\cr
	Version: \tab 2.3.1\cr
	Date: \tab 2017-11-15\cr
	License: \tab GPL (>=2)
}
}
\author{
Peter K Dunn

Maintainer: Peter K Dunn <pdunn2@usc.edu.au>
}
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
	
   Tweedie, M. C. K. (1984).
   An index which distinguishes between some important exponential families.
   \emph{Statistics: Applications and New Directions.
   Proceedings of the Indian Statistical Institute Golden Jubilee International Conference}
   (Eds. J. K. Ghosh and J. Roy), pp. 579--604. Calcutta: Indian Statistical Institute.
}
\keyword{ models }
\examples{
# Generate random numbers
set.seed(987654)
y <- rtweedie( 20, xi = 1.5, mu = 1, phi = 1)
	# With Tweedie index  xi   between 1 and 2, this produces continuous
	# data with exact zeros
x <- rnorm( length(y), 0, 1)  # Unrelated predictor

# With exact zeros, Tweedie index  xi  must be between 1 and 2

# Fit the tweedie distribution; expect  xi  about 1.5
library(statmod)

xi.vec <- seq(1.1, 1.9, by = 0.5)
out <- tweedie.profile( y ~ 1, xi.vec = xi.vec, do.plot = TRUE, verbose = TRUE)

# Fit the glm
require(statmod) # Provides  tweedie  family functions
summary(glm( y ~ x, family = tweedie(var.power = out$xi.max, link.power = 0) ))
}
