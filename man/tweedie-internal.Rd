\name{Tweedie internals}
\alias{dtweedie.dlogfdphi}
\alias{dtweedie.logl}
\alias{dtweedie.logl.saddle}
\alias{dtweedie.logv.bigp}
\alias{dtweedie.logw.smallp}
\alias{dtweedie.interp}
\alias{dtweedie.jw.smallp}
\alias{dtweedie.kv.bigp}
\alias{dtweedie.series.bigp}
\alias{dtweedie.series.smallp}
\alias{stored.grids}
\alias{twcdf}
\alias{twpdf}


\title{Tweedie internal function}

\description{Internal tweedie functions.}

\usage{
	dtweedie.dlogfdphi(y, mu, phi, power)
	dtweedie.logl(phi, y, mu, power)
	dtweedie.logl.saddle( phi, power, y, mu, eps=0)
	dtweedie.logv.bigp( y, phi, power)
	dtweedie.logw.smallp(y, phi, power)
	dtweedie.interp(grid, nx, np, xix.lo, xix.hi,p.lo, p.hi, power, xix)
	dtweedie.jw.smallp(y, phi, power )
	dtweedie.kv.bigp(y, phi, power)
	dtweedie.series.bigp(power, y, mu, phi)
	dtweedie.series.smallp(power, y, mu, phi)
	stored.grids(power)
	twpdf(p, phi, y, mu, exact, verbose, funvalue, exitstatus, relerr, its )
	twcdf(p, phi, y, mu, exact,          funvalue, exitstatus, relerr, its )
}

\arguments{
	\item{y}{the vector of responses}
	\item{power}{the value of \eqn{p}{power} such that the variance is
	\eqn{\mbox{var}[Y]=\phi\mu^p}{var(Y) = phi * mu^power}}
	\item{mu}{the mean}
	\item{phi}{the dispersion}
	\item{grid}{the interpolation grid necessary for the given value of \eqn{p}{power}}
	\item{nx}{the number of interpolation points in the \eqn{\xi}{xi} dimension}
	\item{np}{the number of interpolation points in the \eqn{p}{power} dimension}
	\item{xix.lo}{the lower value of the transformed \eqn{\xi}{xi} value used in the interpolation grid.
	(Note that the value of \eqn{\xi}{xi} is from \eqn{0} to \eqn{\infty}{infty},
	and is transformed such that it is on the range \eqn{0} to \eqn{1}.)}
	\item{xix.hi}{the higher value of the transformed \eqn{\xi}{xi} value used in the interpolation grid.}
	\item{p.lo}{the lower value of \eqn{p} value used in the interpolation grid.}
	\item{p.hi}{the higher value of \eqn{p} value used in the interpolation grid.}
	\item{xix}{the value of the transformed \eqn{\xi}{xi} at which a value is sought.}
	\item{eps}{the offset in computing the variance function in the saddlepoint approximation.
	The default is \code{eps=1/6}
	(as suggested by Nelder and Pregibon, 1987).}
	\item{p}{the Tweedie index parameter}
	\item{exact}{a flag for the FORTRAN to use exact-zeros acceleration algorithmic the calculation (1 means to do so)}
	\item{verbose}{a flag for the FORTRAN: 1 means to be verbose }
	\item{funvalue}{the value of the call returned by the FORTRAN code}
	\item{exitstatus}{the exit status returned by the FORTRAN code}
	\item{relerr}{an estimation of the relative error returned by the FORTRAN code}
	\item{its}{the number of iterations of the algorithm returned by the FORTRAN code}	
}
\details{
These are not to be called by the user.
}

\author{Peter Dunn (\email{pdunn2@usc.edu.au})}
\references{
	Nelder, J. A. and Pregibon, D. (1987).
	An extended quasi-likelihood function
	\emph{Biometrika},
	\bold{74}(2), 221--232.
	doi{10.1093/biomet/74.2.221}
}
\keyword{models}
