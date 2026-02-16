## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----TWexample----------------------------------------------------------------
library(tweedie)
library(statmod)   # For the tweedie() family function, for use in glm()

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
      xlab = expression(italic(y)),
      ylab = "Density function",
      type ="l",
      lwd = 2,
      las = 1)
points(twden[y==0] ~ y[y == 0],
      lwd = 2,
      pch = 19)

## ----TWplotsCDF---------------------------------------------------------------
plot(twdtn ~ y,
     xlab = expression(italic(y)),
     ylab = "Distribution function",
     type = "l",
     las = 1,
     lwd = 2,
     las = 1,
     ylim = c(0, 1) )
points(twdtn[y==0] ~ y[y == 0],
      lwd = 2,
      pch = 19)

## ----TWplotsPDF2--------------------------------------------------------------
tweedie::tweedie_plot(y, xi = xi, mu = mu, phi = phi,
                      ylab = "Density function",
                      xlab = expression(italic(y)),
                      las = 1,
                      lwd = 2)

## ----TWplotsCDF2--------------------------------------------------------------
tweedie::tweedie_plot(y, xi = xi, mu = mu, phi = phi, 
                      ylab = "Distribution function",
                      xlab = expression(italic(y)),
                      las = 1, 
                      lwd = 2, 
                      ylim = c(0, 1), 
                      type = "cdf")

## ----TWqqplot-----------------------------------------------------------------
library(tweedie)

qqnorm( qr <- statmod::qresid(mod.tw),
        las = 1)
qqline(qr,
       col = "grey")

## ----TWprofile----------------------------------------------------------------
out <- tweedie::tweedie_profile(y ~ 1, 
                                xi.vec = seq(1.2, 1.8, by = 0.05), 
                                do.plot = TRUE)

# The estimated power index:
out$xi.max

## ----PSNLoad------------------------------------------------------------------
poison <- read.csv(system.file("extdata", "poison.csv", package = "tweedie"))

# Convert to factors:
poison$Psn <- factor(poison$Psn)
poison$Trmt <- factor(poison$Trmt)

head(poison)

## ----PSNPlotPsn---------------------------------------------------------------
plot(Time ~ Psn,
     data = poison,
     xlab = "Poison type",
     ylab = "Survival time (tens of hours)",
     las = 1)

## ----PSNPlotTmt---------------------------------------------------------------
plot(Time ~ Trmt,
     data = poison,
     xlab = "Treatment type",
     ylab = "Survival time (tens of hours)",
     las = 1)

## ----PSNProfile---------------------------------------------------------------
PSNPrf <- tweedie_profile(Time ~ Trmt + Psn, data = poison, 
                          do.plot = TRUE, xi.vec = seq(2.5, 5.5, length = 11) )

## ----PSNXi--------------------------------------------------------------------
PSNPrf$xi.max

## ----PSNModel-----------------------------------------------------------------
PSNMod <- glm(Time ~ Trmt + Psn, 
              data = poison,
              family=statmod::tweedie(var = 4, link.power = 0))
anova(PSNMod, test = "F")

## ----PSNRes-------------------------------------------------------------------
qr <- statmod::qresid(PSNMod)
qqnorm(qr,
       las = 1)
qqline(qr,
       col = "grey")

