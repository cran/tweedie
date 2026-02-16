test_that("Other functions work OK", {
  # For upstream dependencies
  expect_no_error({
    y <- rgamma(5, 1, 1)
    m1 <- glm(y ~ 1, 
              family = statmod::tweedie(var.power = 2.5, 
                                        link.power = 0))
    logLiktweedie(m1)
  })
  
  
  
  # rtweedie
  expect_no_error({
    rtweedie(3, xi = 1.5, mu = 1:3, phi = c(0.25, 0.5, 1) )
  })
  expect_no_error({
    rtweedie(4, power = 1.5, mu = 1:3, phi = c(0.25, 0.5, 1) )
  })
  
  expect_no_error({
    rtweedie(3, 
             xi = 2.5,
             mu = 1:3,
             phi = c(0.25, 0.5, 1) )
  })
  expect_no_error({
    rtweedie(4, power = 2.5,
             mu = 1:3,
             phi = c(0.25, 0.5, 1) )
  })

  test_that("dtweedie power=1 returns 0 and warns", {
    expect_warning(
      val <- tweedie::dtweedie(1, mu = 2, phi = 2, power = 1),
      regexp = "non-integer x"
    )
    expect_equal(val, 0)
  })
  

  # qtweedie
  expect_no_error({
    qtweedie(c(0.1, 0.3, 0.6), 
             xi = 2.5,
             mu = 1:3,
             phi = c(0.25, 0.5, 1) )
  })
  expect_no_error({
    qtweedie(c(0.1, 0.3, 0.6),
             power = 2.5,
             mu = 1:3,
             phi = c(0.25, 0.5, 1) )
  })
  
  expect_no_error({
    qtweedie(c(0.1, 0.3, 0.6), 
             xi = 1.5,
             mu = 1:3,
             phi = c(0.25, 0.5, 1) )
  })
  expect_no_error({
    qtweedie(c(0.1, 0.3, 0.6),
             power = 1.5,
             mu = 1:3,
             phi = c(0.25, 0.5, 1) )
  })
  

  # tweedie_lambda
    expect_true(
    all( is.na( 
      tweedie_lambda(power = 2.5,
                     mu = 1:3,
                     phi = c(0.25, 0.5, 1) )
    )) )
  

  
  # tweedie_convert
  expect_no_error({
    y <- rgamma(5, 1, 1)
    m1 <- glm(y ~ 1, 
              family = statmod::tweedie(var.power = 1.5, 
                                        link.power = 0))
    tweedie_convert(xi = 1.5,
                    mu = fitted(m1),
                    phi = summary(m1)$dispersion)
  })

  expect_error({
    y <- rgamma(5, 1, 1)
    m1 <- glm(y ~ 1, 
              family = statmod::tweedie(var.power = 2.5, 
                                        link.power = 0))
    tweedie_convert(xi = 2.5,
                    mu = fitted(m1),
                    phi = summary(m1)$dispersion)
  },
  regexp = "xi must be less than 2")

  expect_error({
    y <- rgamma(5, 1, 1)
    m1 <- glm(y ~ 1, 
              family = statmod::tweedie(var.power = 2.5, 
                                        link.power = 0))
    tweedie_convert(p = 2.5,
                    mu = fitted(m1),
                    phi = summary(m1)$dispersion)
  },
  regexp = "power must be less than 2")

  # tweedie_lambda
  expect_no_error({
    tweedie_lambda(power = 1.5,
                    mu = 1:3,
                    phi = c(0.25, 0.5, 1) )
  })

  expect_true(
    all( is.na( 
      tweedie_lambda(power = 2.5,
                   mu = 1:3,
                   phi = c(0.25, 0.5, 1) )
    )) )

})

