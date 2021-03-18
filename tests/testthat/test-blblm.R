fit_r = lm(mpg ~ wt * hp, data = mtcars)
fit = blblm(mpg ~ wt * hp, data = mtcars, m = 3, B = 100)
fit_par = blblm(mpg ~ wt * hp, data = mtcars, m = 3, B = 100, parallel = TRUE)

test_that("blblm object found", {
  expect_s3_class(fit, "blblm")
  expect_s3_class(fit_par, "blblm")
})

test_that("blblm coefficient returns same vector length", {
  coefs_r = coef(fit_r)
  coefs = coef(fit)
  expect_equal(length(coefs), length(coefs_r))
})

test_that("blblm parallelization coefficient returns same vector length", {
  coefs_r = coef(fit_r)
  coefs_par = coef(fit_par)
  expect_equal(length(coefs_par), length(coefs_r))
})

test_that("confidence interval returns correct length", {
  con_int = confint(fit_par, c("wt", "hp"))
  expect_equal(length(con_int), 4)
})

test_that("sigma function works", {
  sigma = sigma(fit_par)
  sigma_ci = sigma(fit_par, confidence = TRUE)
  expect_equal(length(sigma), 1)
  expect_equal(length(sigma_ci), 3)
})

test_that("prediction interval returns correct length", {
  prediction_interval = predict(fit_par, data.frame(wt = c(2.5, 3), hp = c(150, 170)))
  expect_equal(length(prediction_interval), 2)
})
