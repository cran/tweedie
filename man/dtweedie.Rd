\name{Tweedie}
\alias{dtweedie}
\alias{dtweedie.series}
\alias{dtweedie.inversion}
\alias{dtweedie.stable}
\alias{ptweedie}
\alias{ptweedie.series}
\alias{ptweedie.inversion}
\alias{qtweedie}
\alias{rtweedie}
\title{Tweedie Distributions}
\description{Density, distribution function, quantile function and random
     generation for the Tweedie family of distributions}
\usage{
dtweedie(y, power, mu, phi)
dtweedie.series(y, power, mu, phi)
dtweedie.inversion(y, power, mu, phi, exact=TRUE, method)
dtweedie.stable(y, power, mu, phi)
ptweedie(q, power, mu, phi)
ptweedie.series(q, power, mu, phi)
qtweedie(p, power, mu, phi)
rtweedie(n, power, mu, phi)
}
\arguments{
\item{y, q}{vector of quantiles}
\item{p}{vector of probabilities}
\item{n}{the number of observations}
\item{power}{the value of \eqn{p}{power} such that the variance is 
\eqn{\mbox{var}[Y]=\phi\mu^p}{var(Y) = phi * mu^power}}
\item{mu}{the mean}
\item{phi}{the dispersion}
\item{exact}{logical flag; 
if \code{TRUE} (the default),
exact zeros are used with
the \eqn{W}{W}-algorithm of Sidi (1982);
if \code{FALSE},
approximate (asymptotic) zeros are used in place of exact zeros.
Using asymptotic zeros requires less computation but is often less accurate;
using exact zeros can be slower but generally improves accuracy.}
\item{method}{either \code{1}, \code{2} or \code{3},
determining which of three methods to use to compute the density
using the inversion method.
If \code{method} is \code{NULL} (the default),
the optimal method (in terms of relative accuracy) is used, element-by-element of \code{y}.
See the Note in the Details section below}
}
\value{
density (\code{dtweedie}),
probability (\code{ptweedie}),
quantile (\code{qtweedie})
or random sample (\code{rtweedie})
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

    This function evaluates the density or cumulative probability 
    using one of two methods, 
    depending on the combination of parameters. 
    One method is the evaluation of an infinite series.
    The second interpolates some stored values computed from a 
    Fourier inversion technique. 
    
    The function \code{dtweedie.inversion}
    evaluates the density using a Fourier series technique; 
    \code{ptweedie.inversion} does likewise for the cumulative 
    probabilities. 
    The actual code is contained in an external FORTRAN program. 
    Different code is used for \eqn{p>2}{power > 2}
    and for \eqn{1<p<2}{1 < power < 2}. 
    
    The function \code{dtweedie.series} evaluates the density
    using a series expansion; 
    a different series expansion
    is used for \eqn{p>2}{power > 2} and for \eqn{1<p<2}{1 < power < 2}.
    The function \code{ptweedie.series} does likewise for the 
    cumulative probabilities but only for \eqn{1<p<2}{1 < power < 2}. 
    
    The function \code{dtweedie.stable} exploits the link between
    the stable distribution (Nolan, 1997) and Tweedie distributions,
    as discussed in Jorgensen, Chapter 4.
    These are computed using Nolan's algorithm as implemented
    in the  fBasics  package (which is therefore required to use
    the dtweedie.stable function).
    
    The function \code{dtweedie} uses a two-dimensional interpolation procedure to 
    compute the density for some parts of the parameter space from 
    previously computed values found from the series or the 
    inversion. For other parts of the parameter space, 
    the series solution is found. 
    
    \code{ptweedie} returns either the computed series 
    solution or inversion solution.
}
    
\section{Note}{The \code{method}s changed from version 1.4 to 1.5
    (methods 1 and 2 swapped).
    The methods are defined in Dunn and Smyth (2007).
}
\author{Peter Dunn (\email{dunn@usq.edu.au})}
\references{
Dunn, Peter K and Smyth, Gordon K (2007).
Evaluation of Tweedie exponential dispersion model by Fourier inversion
\emph{Statistics and Computing}.
To appear.

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

Nolan, John P (1997).
Numerical calculation of stable densities and distribution functions.
\emph{Communication in Statistics---Stochastic models},
\bold{13}(4). 759--774.

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
\seealso{\code{\link{dtweedie.saddle}} }

\examples{
### Plot a Tweedie density
power <- 2.5
mu <- 1 
phi <- 1 
y <- seq(0, 6, length=500) 
fy <- dtweedie( y=y, power=power, mu=mu, phi=phi) 
plot(y, fy, type="l", lwd=2, ylab="Density")
# Compare to the saddlepoint density
f.saddle <- dtweedie.saddle( y=y, power=power, mu=mu, phi=phi) 
lines( y, f.saddle, col=2 )
legend("topright", col=c(1,2), lwd=c(2,1),
    legend=c("Actual","Saddlepoint") )

### A histogram of Tweedie random numbers
hist( rtweedie( 1000, power=1.2, mu=1, phi=1) )

### An example of the multimodal feature of the Tweedie
### family with power near 1 (from Dunn and Smyth, 2005).
y <- seq(0.001,2,len=1000)
mu <- 1
phi <- 0.1
p <- 1.02
f1 <- dtweedie(y,mu=mu,phi=phi,power=p)
plot(y, f1, type="l", xlab="y", ylab="Density")
p <- 1.05
f2<- dtweedie(y,mu=mu,phi=phi,power=p)
lines(y,f2, col=2)

### Compare series and saddlepoint methods
y <- seq(0.001,2,len=1000)
mu <- 1
phi <- 0.1
p <- 1.02
f.series <- dtweedie.series( y,mu=mu,phi=phi,power=p )
f.saddle <- dtweedie.saddle( y,mu=mu,phi=phi,power=p )

f.all <- c( f.series, f.saddle )
plot( range(f.all) ~ range( y ), xlab="y", ylab="Density", 
  type="n")
lines( f.series ~ y, lty=1, col=1)
lines( f.saddle ~ y, lty=3, col=3)

legend("topright", lty=c(1,3), col=c(1,3),
  legend=c("Series","Saddlepoint") )

} 
\keyword{models}

