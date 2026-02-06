## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(tweedie)

## ----TWexample----------------------------------------------------------------
library(statmod)

set.seed(96)
N <- 25
# Mean of the Poisson (lambda) and Gamma (shape/scale)
lambda <- 1.5
# Generating Compound Poisson-Gamma data manually
y <- replicate(N, {
  n_events <- rpois(1, lambda = lambda)
  if (n_events == 0) 0 else sum(rgamma(n_events, shape = 2, scale = 1))
})

mod.tw <- glm(y ~ 1, 
              family = statmod::tweedie(var.power = 1.5, link.power = 0) )
              # link.power = 0  means the log-link

## ----TWrandom-----------------------------------------------------------------
tweedie::rtweedie(10, xi = 1.1, mu = 2, phi = 1)

## ----TWplotsPDF---------------------------------------------------------------
y <- seq(0, 2, length = 100)
xi <- 1.1
mu <- 0.5
phi <- 0.4

twden <- tweedie::dtweedie(y, xi = xi, mu = mu, phi = phi)  
twdtn <- tweedie::ptweedie(y, xi = xi, mu = mu, phi = phi)

plot( twden[y > 0] ~ y[y > 0], 
      type ="l",
      lwd = 2,
      xlab = expression(italic(y)),
      ylab = "Density function")
points(twden[y==0] ~ y[y == 0],
      lwd = 2,
      pch = 19,
      xlab = expression(italic(y)),
      ylab = "Distribution function")

## ----TWplotsCDF---------------------------------------------------------------
plot(twdtn ~ y,
     type = "l",
      lwd = 2,
     ylim = c(0, 1),
      xlab = expression(italic(y)),
      ylab = "Distribution function")

## ----TWplotsPDF2--------------------------------------------------------------
tweedie::tweedie_plot(y, xi = xi, mu = mu, phi = phi,
                      ylab = "Density function",
                      xlab = expression(italic(y)),
                      lwd = 2)

## ----TWplotsCDF2--------------------------------------------------------------
tweedie::tweedie_plot(y, xi = xi, mu = mu, phi = phi, 
                      ylab = "Distribution function",
                      xlab = expression(italic(y)),
                      lwd = 2, 
                      ylim = c(0, 1), 
                      type = "cdf")

## ----TWqqplot-----------------------------------------------------------------
library(tweedie)

qqnorm( statmod::qresid(mod.tw) )

## ----TWprofile----------------------------------------------------------------
out <- tweedie::tweedie.profile(y ~ 1, 
                                xi.vec = seq(1.2, 1.8, by = 0.05), 
                                do.plot = TRUE)

# The estimated power index:
out$xi.max

