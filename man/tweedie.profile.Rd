\name{tweedie.profile}
\alias{tweedie.profile}
\title{Tweedie Distributions: mle estimation of p}
\description{Maximum likelihood estimation of the Tweedie index parameter \eqn{p}{power}.}
\usage{tweedie.profile(formula, p.vec=NULL, link.power=0, data, weights, offset, fit.glm=FALSE,
do.smooth=TRUE, do.plot=FALSE, do.ci=do.smooth,
eps=1/6, do.points=do.plot, method="inversion", conf.level=0.95, 
phi.method=ifelse(method == "saddlepoint", "saddlepoint", "mle"), verbose=FALSE)}
\arguments{
   \item{formula}{a formula expression as for other regression models and generalized linear models, 
   of the form \code{response ~ predictors}. 
   For details, 
   see the documentation for \code{\link{lm}}, 
   \code{\link{glm}} and \code{\link{formula}}}
   \item{p.vec}{a vector of \code{p} values for consideration.
   The values must all be larger than one
   (if the response variable has exact zeros,
   the values must all be between one and two).
   If \code{NULL} (the default),
   \code{p.vec} is set to 
   \code{seq(1.2, 1.8, by=0.1)} if the 
   response contains any zeros,
   or
   \code{seq(1.5, 5, by=0.25)} if the 
   response contains no zeros.
   See the DETAILS section below for further details.}
   \item{link.power}{the power link function to use.
   These link functions \eqn{g(\cdot)}{g()} are of the form 
   \eqn{g(\eta)=\eta^{\rm link.power}}{g(eta) = eta^link.power},
   and the special case of \code{link.power=0} (the default)
   refers to the logarithm link function.
   See the documentation for 
   \code{\link[statmod]{tweedie}} also.}
   \item{data}{an optional data frame, list or environment 
  	(or object coercible by \code{as.data.frame} to a data frame) 
  	containing the variables in the model.  
  	If not found in \code{data}, 
  	the variables are taken from \code{environment(formula)}, 
  	typically the environment from which \code{glm} is called.}
  	\item{weights}{an optional vector of weights to be used in the fitting
  	process.  Should be \code{NULL} or a numeric vector.}
  	\item{offset}{this can be used to specify an \emph{a priori}
  	 known component to be included in the linear predictor during fitting.  
  	 This should be \code{NULL} or a numeric vector of length either one or
  	 equal to the number of cases. 
  	 One or more \code{offset} terms can
  	 be included in the formula instead or as well, 
  	 and if both are specified their sum is used.  
  	 See \code{\link{model.offset}}.}
   \item{fit.glm}{logical flag.
   If \code{TRUE},
   the Tweedie generalized linear model is fitted using the value of \eqn{p}{p}
   found by the profiling function.
   If \code{FALSE} (the default),
   no model is fitted.}
   \item{do.smooth}{logical flag.
   If \code{TRUE} (the default),
   a spline is fitted to the data to smooth the profile likelihood plot.
   If \code{FALSE},
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
   no confidence interval is computed.
   By default,
   \code{do.ci} is the same value as \code{do.smooth},
   since a confidence interval will only be accurate if
   smoothing has been performed.
   Indeed,
   if \code{do.smooth=FALSE},
   confidence intervals are never computed and
   \code{do.ci} is forced to \code{FALSE} if it is given as \code{TRUE}.}
   \item{eps}{the offset in computing the variance function.
   The default is \code{eps=1/6}
   (as suggested by Nelder and Pregibon, 1987).
   Note \code{eps} is ignored unless the 
   \code{method="saddlepoint"}
   as it makes no sense otherwise.}
   \item{do.points}{plot the points on the plot where the
   (log-) likelihood is computed for the given values of \code{p};
   defaults to the same value as \code{do.plot}}
   \item{method}{the method for computing the (log-) likelihood.
   One of
   \code{"series"},
   \code{"inversion"} (the default),
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
   \item{verbose}{the amount of feedback requested:
   \code{0} or \code{FALSE} means minimal feedback (the default), 
   \code{1} or \code{TRUE} means some feedback,
   or \code{2} means to show all feedback.
   Since the function can be slow and sometimes problematic,
   feedback can be good;
   but it can also be unnecessary when one knows all is well.}
}
\value{
   The main purpose of the function is to estimate the value
   of the Tweedie index parameter, \eqn{p}{p},
   which is produced by the output list as \code{p.max}.
   Optionally (if \code{do.plot=TRUE}),
   a plot is produced that shows the profile log-likelihood
   computed at each value in \code{p.vec}
   (smoothed if \code{do.smooth=TRUE}).
   This function can be tempermental 
   (for theoretical reasons involved in numerically computing the density),
   and this plot shows the values of \eqn{p}{p} requested on the
   horizontal axis (using \code{\link{rug}}); 
   there may be fewer points on the plot,
   since the likelihood some values of \eqn{p}{p} requested
   may have returned \code{NaN}, \code{Inf} or \code{NA}.
   
   A list containing the components:
   \code{y} and \code{x}
   (such that \code{plot(x,y)} (partially)
   recreates the profile likelihood plot);
   \code{ht} (the height of the nominal confidence interval);
   \code{L} (the estimate of the (log-) likelihood at each given value of \code{p});
   \code{p} (the \code{p}-values used);
   \code{phi} (the computed values of \code{phi} at the values in \code{p});
   \code{p.max} (the estimate of the mle of \code{p});
   \code{L.max} (the estimate of the (log-) likelihood at \code{p.max});
   \code{phi.max} (the estimate of \code{phi} at \code{p.max});
   \code{ci} (the lower and upper limits of the confidence interval for \code{p});
   \code{method} (the method used for estimation: \code{series}, \code{inversion}, 
   \code{interpolation} or \code{saddlepoint});
   \code{phi.method} (the method used for estimation of \code{phi}:
   \code{saddlepoint} or \code{phi}).
   
   If \code{glm.fit} is \code{TRUE},
   the list also contains a component \code{glm.obj},
   a \code{glm} object  for the fitted Tweedie generalized linear model.
}
\note{
   The estimates of \code{p}
   and \code{phi} are printed.
   The result is printed invisibly.
   
   If the response variable has any exact zeros,
   the values in \code{p.vec}
   must all be between one and two.
   
   The function is sometimes unstable and may fail.
   It may also be very slow.
   One solution is to change the method.
   The default is \code{method="inversion"} (the default);
   then try \code{method="series"},
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
   \code{do.smooth=FALSE},
   \code{do.ci=FALSE}.
   If this is successful,
   a larger vector \code{p.vec}
   and smoothing can be used.
}
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
\author{
   Peter Dunn (\email{pdunn2@usc.edu.au})
}
\references{
   Dunn, P. K. \& Smyth, G. K. (2008).
   Evaluation of Tweedie exponential dispersion model densities by Fourier inversion.
   \emph{Statistics and Computing}, 
   \bold{18}, 73--86.
   
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
   \code{\link[statmod]{tweedie}}
}

\examples{
library(statmod) # Needed to use  tweedie.profile
# Generate some fictitious data
test.data <- rgamma(n=200, scale=1, shape=1)
# The gamma is a Tweedie distribution with power=2;
# let's see if p=2 is suggested by  tweedie.profile:
\dontrun{
	out <- tweedie.profile( test.data ~ 1, 
		p.vec=seq(1.5, 2.5, by=0.2) )
	out$p.max
	out$ci
}
}

\keyword{models}
