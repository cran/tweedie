\name{tweedie.profile}
\alias{tweedie.profile}
\title{Tweedie Distributions: mle estimation of p}
\description{Maximum likelihood estimation of the Tweedie index parameter \eqn{p}{power}.}
\usage{tweedie.profile(formula, p.vec, smooth=FALSE, do.plot=FALSE, do.ci=smooth,
eps=1/6, method="series", conf.level=0.95, 
phi.method=ifelse(method == "saddlepoint", "saddlepoint", "mle"))}
\arguments{
\item{formula}{a formula expression as for other regression models and generalized linear models, 
of the form \code{response \eqn{\sim}{~} predictors}. 
For details, 
see the documentation for \code{\link{lm}}, 
\code{\link{glm}} and \code{\link{formula}}}
\item{p.vec}{a vector of \code{p} values for consideration.
The values must all be larger than one
(if the response variable has exact zeros,
the values must all be between one and two).
See the DETAILS section below for further details.}
\item{smooth}{logical flag.
If \code{TRUE},
a spline is fitted to the data to smooth the profile likelihood plot.
If \code{FALSE} (the default),
no smoothing is used 
(and the function is quicker).
\bold{Note} that \code{p.vec} must contain \emph{at least five points}
for smoothing to be allowed.}
\item{do.plot}{logical flag.
If \code{TRUE},
a plot of the profile likelihood is produce.
If \code{FALSE} (the default),
no plot is produced.}
\item{do.ci}{logical flag.
If \code{TRUE},
the nominal 100*\code{conf.level}
is computed.
If \code{FALSE},
the confidence interval is not computed.
By default,
\code{do.ci} is the same value as \code{smooth},
since a confidence interval will only be accurate if
smoothing has been performed.}
\item{eps}{the offset in computing the variance function.
The default is \code{eps=1/6}
(as suggested by Nelder and Pregibon, 1987).
Note that \code{eps} is ignored unless the 
\code{method="saddlepoint"}
as it makes no sense.}
\item{method}{the method for computing the (log-) likelihood.
One of
\code{"series"} (the default),
\code{"inversion"},
\code{"interpolation"}
or
\code{"saddlepoint"}.
If there are any troubles using this function,
often a change of method will fix the problem.
Note that \code{method="saddlepoint"}
is only an approximate method for computing the (log-) likelihood.
Using \code{method="interpolation"}
may produce a jump in the profile likelihood as it changes computational regimes.}
\item{conf.level}{the confidence level for the computation of the nominal
confidence interval.
The default is \code{conf.level=0.95}.}
\item{phi.method}{the method for estimating \code{phi},
one of
\code{"saddlepoint"}
or
\code{"mle"}.
A maximum likelihood estimate is used unless
\code{method="saddlepoint"},
when the saddlepoint approximation method is used.
Note that using 
\code{phi.method="saddlepoint"}
is equivalent to using the mean deviance estimator of \code{phi}.
}
}
\value{
A list containing the components:
\code{y} and \code{x}
(such that \code{plot(x,y)} (partially)
recreates the profile likelihood plot);
\code{ht} (the height of the nominal confidence interval);
\code{L} (the estimate of the (log-) likelihood at each given value of \code{p});
\code{p} (the \code{p}-values used);
\code{p.max} (the estimate of the mle of \code{p});
\code{L.max} (the estimate of the (log-) likelihood at \code{p.max});
\code{phi} (the estimate of \code{phi} at \code{p.max});
\code{ci} (the lower and upper limits of the confidence interval for \code{p});
\code{method} (the method used).
}
\note{The estimates of \code{p}
and \code{phi} are printed.
The result is printed invisibly.

If the response variable has any exact zeros,
the values in \code{p.vec}
must all be between one and two.

The function is sometimes unstable and may fail.
It may also be very slow.
One solution is to change the method.
The default is \code{method="series"} (the default);
then try \code{method="inversion"},
\code{method="interpolation"}
and
\code{method="saddlepoint"}
in that order.
Note that 
\code{method="saddlepoint"}
is an approximate method only.
Also make sure the values in \code{p.vec}
are suitable for the data  
(see above paragraph).

It is recommended that for the first use with a data set,
use \code{p.vec} with only a small number of values
and set
\code{smooth=FALSE},
\code{do.ci=FALSE}.
If this is successful,
a larger vector \code{p.vec}
and smoothing can be used.}
\details{
For each value in \code{p.vec},
the function computes an estimate of \code{phi}
and then computes the value of the log-likelihood for these parameters.
The plot of the log-likelihood against \code{p.vec} 
allows the maximum likelihood value of \code{p}
to be found.
Once the value of \code{p} is found,
the distribution within the class of Tweedie distribution is identified.
}
\author{Peter Dunn (\email{dunn@usq.edu.au})}
\references{
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

Nelder, J. A. and Pregibon, D. (1987).
An extended quasi-likelihood function.
\emph{Biometrika}
\bold{74}(2),
221--232.

Tweedie, M. C. K. (1984).
An index which distinguishes between some important exponential families.
\emph{Statistics: Applications and New Directions.
Proceedings of the Indian Statistical Institute Golden Jubilee International Conference}
(Eds. J. K. Ghosh and J. Roy), pp. 579-604. Calcutta: Indian Statistical Institute.
}
\seealso{
\code{\link{dtweedie}}
\code{\link{dtweedie.saddle}}
}
\examples{
library(statmod) # Needed to use  tweedie.profile
# Generate some fictitious data
test.data <- rgamma(n=200, scale=1, shape=1)
# The gamma is a Tweedie distribution with power=2;
# let's see if the profile plot shows this
out <- tweedie.profile( test.data ~1, p.vec=seq(1.7, 2.3, length=10),
       do.plot=TRUE, method="interpolation", smooth=TRUE, do.ci=TRUE)
# And plot the points to see how the smooth went
points( out$p, out$L)
}

\keyword{models}
