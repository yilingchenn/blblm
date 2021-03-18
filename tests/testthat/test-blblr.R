mydata = read.csv("https://stats.idre.ucla.edu/stat/data/binary.csv")
fit_r = glm(admit ~ gre + gpa + rank, data = mydata, family = "binomial")
fit = blblr(admit ~ gre + gpa + rank, data = mydata, m = 3, B = 100)
# fit_par = blblr(admit ~ gre + gpa + rank, data = mydata, m = 3, B = 100, parallel = TRUE)
newdata = with(mydata, data.frame(gre = mean(gre), gpa = mean(gpa), rank = factor(1:4)))

test_that("blblr object found", {
  expect_s3_class(fit, "blblr")
  # expect_s3_class(fit_par, "blblr")
})

test_that("blblr coefficient returns same vector length", {
  coefs_r = coef(fit_r)
  coefs = coef(fit)
  expect_equal(length(coefs), length(coefs_r))
})

test_that("confidence interval returns correct length", {
  con_int = confint(fit, c("gre", "gpa", "rank"))
  con_int2 = confint(fit_r, c("gre", "gpa", "rank"))
  expect_equal(length(con_int), length(con_int2))
})

test_that("prediction interval returns correct length", {
  prediction_interval = predict(fit, data.frame(gre = c(580, 588), gpa = c(3.13, 3.67), rank = c(1.00, 2.48)))
  expect_equal(length(prediction_interval), 2)
})