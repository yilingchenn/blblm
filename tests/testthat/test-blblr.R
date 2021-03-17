fit_r = glm(mpg ~ wt*hp, data = mtcars)
coefs_r = coef(fit_r)
fit = blblr(mpg ~ wt*hp, data = mtcars)
fit_par = blblr(mpg ~ wt*hp, data = mtcars, parallel = TRUE)
coefs = coef(fit)
coefs_par = coef(fit_par)

test_that("blblr coefficient returns same vector length", {
  expect_s3_class(fit, "blblr")
  expect_equal(length(coefs), length(coefs_r))
})

test_that("blblr parallelization coefficient returns same vector length", {
  expect_s3_class(fit, "blblr")
  expect_equal(length(coefs_par), length(coefs_r))
})

test_that("return blblr object", {
  expect_equal(length(fit), length(fit_r))
})

test_that("confidence interval", {

})