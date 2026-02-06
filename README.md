
# tweedie <a href="https://CRAN.R-project.org/package=tweedie/"><img src="man/figures/tweedie_hex.png" align="right" height="138" /></a>

<!-- README.md is generated from README.Rmd. Please edit that file -->

<!-- badges: start -->

<!-- ![R-CMD-check](https://github.com/rstudio/rmarkdown/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/rstudio/rmarkdown/actions/workflows/R-CMD-check.yaml)
[![CRAN release](https://www.r-pkg.org/badges/version/rmarkdown)](https://cran.r-project.org/package=rmarkdown)
[![Codecov test coverage](https://codecov.io/gh/rstudio/rmarkdown/branch/main/graph/badge.svg)](https://app.codecov.io/gh/rstudio/rmarkdown?branch=main) -->

<!-- badges: end -->

The **tweedie** package allows likelihood computations for Tweedie
distributions.

Apart from special cases (the normal, Poisson, gamma, inverse Gaussian
distributions), Tweedie distributions do not have closed-form density
functions or distribution functions. This package uses fast numerical
algorithms (infinite oscillation integrals; infinite series) to evaluate
the Tweedie density functions and distribution functions.

## Installation

You can install the development version of **tweedie** from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("PeterKDunn/tweedie")
```

## Tweedie distributions

Tweedie distributions are exponential dispersion models, with a mean
$\mu$ and a variance $\phi \mu^\xi$, for some dispersion parameter
$\phi > 0$ and a power index $\xi$ (sometimes called $p$) that uniquely
defines the distribution within the Tweedie family (for all real values
of $\xi$ **not** between 0 and 1).

Special cases of the Tweedie distributions are:

- the *normal* distribution, with $\xi = 0$ (i.e., the variance is
  $\phi$ and not related to the mean);
- the *Poisson* distribution, with $\xi = 1$ and $\phi = 1$ (i.e., the
  variance is the same as the mean);
- the *gamma* distribution, with $\xi = 2$; and
- the *inverse Gaussian* distribution, with $\xi = 3$.

For all other values of $\xi$, the probability functions and
distribution functions have no closed forms.

For $\xi < 1$, applications are limited (non-existent so far?), but have
support on the entire real line and $\mu > 0$.

For $1 < \xi < 2$, Tweedie distributions can be represented as a Poisson
sum of gamma distributions. These distributions are continuous for
$Y > 0$ but have a discrete mass at $Y = 0$.

For $\xi \ge 2$, the distributions have support on the positive reals.

The vignette contains examples.

<!-- REMEMBER to commit and push the resulting figure files, so they display on GitHub and CRAN. -->
