
test_that("Deprecated functions work", {
  expect_warning(
    tweedie.dev(y=1, mu = 1, power = 6)  
  )
  expect_warning(
    ptweedie.inversion(1, mu = 1, phi = 1, power = 6, verbose=FALSE)  
  )
  expect_warning(
    dtweedie.inversion(1, mu = 1, phi = 1, power = 6, verbose=FALSE)  
  )
  expect_warning(
    dtweedie.series(1, mu = 1, phi = 1, power = 6)  
  )
  expect_warning(
    ptweedie.series(1, mu = 1, phi = 1, power = 1.5)  
  )
  expect_warning(
    tweedie.plot(1, mu = 1, phi = 1, power = 1.5)  
  )
  expect_warning({
    set.seed(1000)
    y <- rgamma(30, shape = 1.2, rate = 1.1)
    tweedie.profile(formula = (y ~ 1), 
                    data = data.frame(y = y ),
                    p.vec = seq(1.5, 2.5, length = 10),
                    link.power = 0, # Log link
                    do.smooth = TRUE, 
                    do.plot = FALSE, 
                    do.ci = FALSE, 
                    method = "interpolation",
                    do.points = FALSE,
                    phi.method = "mle")},
    regexp = "deprecated"
  )
  expect_warning(
    tweedie.plot(y = seq(0, 10, by = 1),
                 mu = 1,
                 phi = 1,
                 power = 1.5)
  )
  expect_warning({
    set.seed(1000)
    y <- rgamma(10, shape = 1.4, scale = 1.5)
    mod1 <- glm(y ~ 1, family=statmod::tweedie(link.power = 0, 
                                               var.power = 2) )
    AICtweedie( mod1, 
                dispersion = NULL,
                k = 2,
                verbose = FALSE)
  }
  )
})  
