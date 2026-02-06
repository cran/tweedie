test_that("Same values as non-central chi-sq (A)", {
  mu <- 1
  phi <- 4
  p <- 1.5
  y <- c(0.001, 0.01, 0.05, 0.1, 0.5, 0.75,
         seq(01, 10, by = 1),
         15, 20, 50, 100) 
  expect_equal(
      ptweedie(y, mu = mu, phi = phi, power = p),
      pchisq(y, df = 0, ncp = 1)
  )

  expect_equal(
    dtweedie(y, mu = mu, phi = phi, power = p),
    dchisq(y, df = 0, ncp = 1)
  )

  # For the inversion, avoid very small y values  
  y <- c(         0.01, 0.05, 0.1, 0.5, 0.75,
         seq(01, 10, by = 1),
         15, 20, 50, 100) 
  
  expect_equal(
    dtweedie_inversion(y, mu = mu, phi = phi, power = p),
    dchisq(y, df = 0, ncp = 1)
  )
  expect_equal(
    ptweedie_inversion(y, mu = mu, phi = phi, power = p),
    pchisq(y, df = 0, ncp = 1)
  )

})



test_that("Same values as inverse Gaussian", {
  library(statmod)
  mu <- 1.4
  phi <- 0.74
  p <- 3
  y <- c(0.001, 0.01, 0.05, 0.1, 0.5, 0.75,
         seq(01, 10, by = 1),
         15, 20, 50, 100) 

  expect_equal(
    statmod::dinvgauss(y, mean = mu, dispersion = phi), 
    dtweedie(y, mu = mu, phi = phi, power = p)
  )
  
  expect_equal(
    statmod::dinvgauss(y, mean = mu, dispersion = phi), 
    dtweedie_inversion(y, mu = mu, phi = phi, power = p, IGexact = FALSE )
  )

  expect_equal(
    statmod::pinvgauss(y, mean = mu, dispersion = phi), 
    ptweedie(y, mu = mu, phi = phi, power = p)
  )

  expect_equal(
    statmod::pinvgauss(y, mean = mu, dispersion = phi), 
    ptweedie_inversion(y, mu = mu, phi = phi, power = p, IGexact = FALSE )
  )

  
  mu <- 0.4
  phi <- 2
  p <- 3
  y <- c(0.001, 0.01, 0.05, 0.1, 0.5, 0.75,
         seq(01, 10, by = 1),
         15, 20, 50, 100) 
  
  expect_equal(
    statmod::dinvgauss(y, mean = mu, dispersion = phi), 
    dtweedie(y, mu = mu, phi = phi, power = p)
  )
  
  expect_equal(
    statmod::dinvgauss(y, mean = mu, dispersion = phi), 
    dtweedie_inversion(y, mu = mu, phi = phi, power = p, IGexact = FALSE )
  )
  
  expect_equal(
    statmod::pinvgauss(y, mean = mu, dispersion = phi), 
    ptweedie(y, mu = mu, phi = phi, power = p)
  )
  
  expect_equal(
    statmod::pinvgauss(y, mean = mu, dispersion = phi), 
    ptweedie_inversion(y, mu = mu, phi = phi, power = p, IGexact = FALSE )
  )
  
  
  
  mu <- 4
  phi <- 2
  p <- 3
  y <- c(0.001, 0.01, 0.05, 0.1, 0.5, 0.75,
         seq(01, 10, by = 1),
         15, 20, 50, 100) 
  
  expect_equal(
    statmod::dinvgauss(y, mean = mu, dispersion = phi), 
    dtweedie(y, mu = mu, phi = phi, power = p)
  )
  
  expect_equal(
    statmod::dinvgauss(y, mean = mu, dispersion = phi), 
    dtweedie_inversion(y, mu = mu, phi = phi, power = p, IGexact = FALSE )
  )
  
  expect_equal(
    statmod::pinvgauss(y, mean = mu, dispersion = phi), 
    ptweedie(y, mu = mu, phi = phi, power = p)
  )
  
  expect_equal(
    statmod::pinvgauss(y, mean = mu, dispersion = phi), 
    ptweedie_inversion(y, mu = mu, phi = phi, power = p, IGexact = FALSE )
  )  
})




