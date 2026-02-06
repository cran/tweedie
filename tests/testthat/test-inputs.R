test_that("Negative y-values return 0 for density", {
  expect_equal(
      dtweedie(-1, mu = 1, phi = 1, power = 3.5),
      0
  )
})


test_that("Negative y-values return 0 for DF", {
  fy <- ptweedie(-1, mu = 1, phi = 1, power = 3.5)
  expect_equal(fy, 0)

  fy <- ptweedie(-1, mu = 1, phi = 1, power = 1.5)
  expect_equal(fy, 0)
})



test_that("Unequal-size input vectors", {
  
  expect_error(  
    ptweedie(c(0, 1, 4), mu = c(1, 3), phi = 1, power = 3.5),
    "mu must be scalar, or the same length as the first input."
  )
  
  expect_error(  
    ptweedie(c(0, 1, 4), mu = 1, phi = c(1, 2), power = 3.5),
    "phi must be scalar, or the same length as the first input."
  )
  
  expect_error(
    ptweedie(c(0, 1, 4), mu = 1, phi = 2, power = c(1.5, 3.5) ),
    "The Tweedie index parameter \\(power\\) must be a single value."
  )
  
  expect_message(
    ptweedie(1, mu = 1, phi = 1, power = 1.5, xi = 1.2),
    "Both  xi  and  power  are given; using the given value of  xi"
  )

  expect_error(
    ptweedie(c(0, 1, 4), mu = 1, phi = 2, xi = c(1.5, 3.5) ),
    "The Tweedie index parameter \\(xi\\) must be a single value."
  )
  
  expect_error(  
    dtweedie(c(0, 1, 4), mu = c(1, 3), phi = 1, power = 3.5),
    "mu must be scalar, or the same length as the first input."
  )
  
  expect_error(  
    dtweedie(c(0, 1, 4), mu = 1, phi = c(1, 2), power = 3.5),
    "phi must be scalar, or the same length as the first input."
  )
  
  expect_error(
    dtweedie(c(0, 1, 4), mu = 1, phi = 2, power = c(1.5, 3.5) ),
    "The Tweedie index parameter \\(power\\) must be a single value."
  )
  expect_error(
    dtweedie(c(0, 1, 4), mu = 1, phi = 2, xi = c(1.5, 3.5) ),
    "The Tweedie index parameter \\(xi\\) must be a single value."
  )
  
  expect_message(
    dtweedie(1, mu = 1, phi = 1, power = 1.5, xi = 1.2),
    "Both  xi  and  power  are given; using the given value of  xi"
  )
})


test_that("Invalid parameters cause error in dtweedie", {
  # 1. Test for bad mu value (mu <= 0)
  # The error message is the string "mu must be positive."
  expect_error(
    dtweedie(1, mu = -1, phi = 1, power = 3.5), 
    regexp = "mu must be positive." 
  )
  
  # 2. Test for bad phi value (phi <= 0)
  # Assuming the error message is "phi must be positive."
  expect_error(
    dtweedie(1, mu = 1, phi = -1, power = 1.5),
    regexp = "phi must be positive."
  )
})




test_that("Invalid parameters cause error in ptweedie", {
  
  # 1. Test for bad mu value (mu <= 0)
  # The error message is the string "mu must be positive."
  expect_error(
    ptweedie(1, mu = -1, phi = 1, power = 3.5), 
    regexp = "mu must be positive." 
  )
  
  # 2. Test for bad phi value (phi <= 0)
  # Assuming the error message is "phi must be positive."
  expect_error(
    ptweedie(1, mu = 1, phi = -1, power = 1.5),
    regexp = "phi must be positive."
  )
})


test_that("Invalid parameters cause error in rtweedie", {
  
  # 1. Test for bad mu value (mu <= 0)
  # The error message is the string "mu must be positive."
  expect_error({
    set.seed(1000)
    rtweedie(2, mu = -1, phi = 1, power = 3.5)}, 
    regexp = "mu must be positive." 
  )
  
  # 2. Test for bad phi value (phi <= 0)
  # Assuming the error message is "phi must be positive."
  expect_error({
    set.seed(1000)
    rtweedie(1, mu = 1, phi = -1, power = 1.5)},
    regexp = "phi must be positive."
  )

  # 3. Test for bad n value (n > 0)
  expect_error({
    set.seed(1000)
    rtweedie(0, mu = 1, phi = -1, power = 1.5)},
    regexp = "n must be a positive integer."
  )
})


test_that("Invalid parameters cause error in qtweedie", {
  
  # 1. Test for bad mu value (mu <= 0)
  expect_error(
    qtweedie(1, mu = -1, phi = 1, power = 1.5), 
    regexp = "mu must be positive." 
  )
  expect_error(
    qtweedie(1, mu = -1, phi = 1, power = 3.5), 
    regexp = "mu must be positive." 
  )
  
  # 2. Test for bad phi value (phi <= 0)
  expect_error(
    qtweedie(1, mu = 1, phi = -1, power = 1.5),
    regexp = "phi must be positive."
  )
  expect_error(
    qtweedie(1, mu = 1, phi = -1, power = 3.5),
    regexp = "phi must be positive."
  )
  
  # 3. Test for q less than zero
  expect_error(
    qtweedie(-1, mu = 1, phi = -1, power = 1.5),
    regexp = "quantiles must be between 0 and 1."
  )
  expect_error(
    qtweedie(-1, mu = 1, phi = -1, power = 3.5),
    regexp = "quantiles must be between 0 and 1."
  )
  
  
  # 4. Test for q greater than zero
  expect_error(
    qtweedie(2, mu = 1, phi = -1, power = 1.5),
    regexp = "quantiles must be between 0 and 1."
  )
  expect_error(
    qtweedie(2, mu = 1, phi = -1, power = 3.5),
    regexp = "quantiles must be between 0 and 1."
  )
  
  
  # 5. Test for q greater than zero, less than 0  
  expect_error(
    qtweedie(c(-1, 2), mu = 1, phi = -1, power = 1.5),
    regexp = "quantiles must be between 0 and 1."
  )
  expect_error(
    qtweedie(c(-1, 2), mu = 1, phi = -1, power = 3.5),
    regexp = "quantiles must be between 0 and 1."
  )
})


