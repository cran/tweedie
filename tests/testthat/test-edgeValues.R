test_that("Zero first inputs returns correct values", {
  
  # 1. Test dtweedie(0, ...)
  expect_equal(
    dtweedie(0, mu = 1, phi = 1, power = 1.5), 
    exp( -1^(2-1.5)/((1 * (2 - 1.5))))
  )
  
  expect_equal(
    dtweedie(0, mu = 1, phi = 1, power = 3.5), 
    0
  )

  # 2. Test ptweedie(0, ...)
  expect_equal(
    ptweedie(0, mu = 1, phi = 1, power = 1.5), 
    exp( -1^(2-1.5)/((1 * (2 - 1.5))))
  )

  expect_equal(
    ptweedie(0, mu = 1, phi = 1, power = 3.5), 
    0
  )
})





