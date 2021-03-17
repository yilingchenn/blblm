fit_r = lm(mpg ~ wt*hp, data = mtcars)
coefs_r = coef(fit_r)
sigma_r = sigma(fit_r)
fit = blblm(mpg ~ wt*hp, data = mtcars)
fit_par = blblm(mpg ~ wt*hp, data = mtcars, parallel = TRUE)
coefs = coef(fit)
sigma = sigma(fit)
coefs_par = coef(fit_par)

test_that("blblm coefficient returns same vector length", {
  expect_s3_class(fit, "blblm")
  expect_equal(length(coefs), length(coefs_r))
})

test_that("blblm parallelization coefficient returns same vector length", {
  expect_s3_class(fit, "blblm")
  expect_equal(length(coefs_par), length(coefs_r))
})

test_that("sigma function works", {
  expect_s3_class(fit, "blblm")
  expect_equal(length(sigma), length(sigma_r))
})

test_that("return blblm object", {
  expect_equal(length(fit), length(fit2))
})

test_that("confidence interval", {

})