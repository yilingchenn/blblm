#' @import purrr
#' @import furrr
#' @import stats
#' @import future
#' @import tidyverse
#' @import utils
#' @importFrom magrittr %>%
#' @details
#' Linear Regression with Little Bag of Bootstraps
"_PACKAGE"

## quiets concerns of R CMD check re: the .'s that appear in pipelines
# from https://github.com/jennybc/googlesheets/blob/master/R/googlesheets.R
utils::globalVariables(c("."))

#' @param formula linear regression formula
#' @param data dataset
#' @param m number of split, default = 10
#' @param B number of bootstrap, default = 5000
#' @param parallel estimation being calculated parallely, default = FALSE
#' @param threads number of workers using to run the estimation, default = 4
#' @name blblm
#' @title blblm
#' @export
blblm <- function(formula, data, m = 10, B = 5000, parallel = FALSE, threads = 4) {
  data_list <- split_data(data, m)
  if (parallel){
    suppressWarnings(plan(multiprocess, workers = threads))
    options(future.rng.onMisuse = "ignore")
    estimates <- future_map(
      data_list,
      ~ lm_each_subsample(formula = formula, data = ., n = nrow(data), B = B))
  }else{
    estimates <- map(
      data_list,
      ~ lm_each_subsample(formula = formula, data = ., n = nrow(data), B = B))
  }
  res <- list(estimates = estimates, formula = formula)
  class(res) <- "blblm"
  invisible(res)
}


#' @param data dataset
#' @param m number of splits
#' split data into m parts of approximated equal sizes
split_data <- function(data, m) {
  idx <- sample.int(m, nrow(data), replace = TRUE)
  data %>% split(idx)
}


#' @param formula linear regression formula
#' @param data dataset
#' @param n number of rows in each subset
#' @param B number of bootstraps
#' compute the estimates
lm_each_subsample <- function(formula, data, n, B) {
  # drop the original closure of formula,
  # otherwise the formula will pick a wrong variable from the global scope.
  environment(formula) <- environment()
  m <- model.frame(formula, data)
  X <- model.matrix(formula, m)
  y <- model.response(m)
  # nrow(data) = n/m, the split data
  replicate(B, lm1(X, y, n), simplify = FALSE)
}

#' @param X design matrix of dimension n * p
#' @param y vector of observations of length n, or a matrix with n rows
#' @param n number of rows in each subset
#' compute the regression estimates for a blb dataset
lm1 <- function(X, y, n) {
  freqs <- as.vector(rmultinom(1, n, rep(1, nrow(X))))
  fit <- lm.wfit(X, y, freqs)
  list(coef = blbcoef(fit), sigma = blbsigma(fit))
}


#' @param fit linear regression model
#' compute the coefficients from fit
blbcoef <- function(fit) {
  coef(fit)
}


#' @param fit linear regression model
#' compute sigma from fit
blbsigma <- function(fit) {
  p <- fit$rank
  e <- fit$residuals
  w <- fit$weights
  sqrt(sum(w * (e^2)) / (sum(w) - p))
}


#' @param x dataset for blblr model
#' @export
#' @method print blblm
print.blblm <- function(x, ...) {
  cat("blblm model:", capture.output(x$formula))
  cat("\n")
}


#' @param object blblm object
#' @param confidence confidence interval, default = FALSE
#' @param level level of confidence interval, default = 0.95
#' @param parallel estimation being calculated parallely, default = FALSE
#' @param threads number of workers using to run the estimation, default = 4
#' @return variance of blblm estimates
#' @name sigma.blblm
#' @export
#' @method sigma blblm
sigma.blblm <- function(object, confidence = FALSE, level = 0.95, parallel = FALSE, threads = 4, ...) {
  est <- object$estimates
  sigma <- mean(map_dbl(est, ~ mean(map_dbl(., "sigma"))))
  if (confidence) {
    alpha <- 1 - 0.95
    limits <- est %>%
      map_mean(~ quantile(map_dbl(., "sigma"), c(alpha / 2, 1 - alpha / 2)), parallel, threads) %>%
      set_names(NULL)
    return(c(sigma = sigma, lwr = limits[1], upr = limits[2]))
  } else {
    return(sigma)
  }
}

#' @param object blblm model
#' @param parallel estimation being calculated parallely, default = FALSE
#' @param threads number of workers using to run the estimation, default = 4
#' @export
#' @method coef blblm
coef.blblm <- function(object, parallel = FALSE, threads = 4,...) {
  est <- object$estimates
  map_mean(est, ~ map_cbind(., "coef") %>% rowMeans(), parallel, threads)
}

#' @param object blblm object
#' @param parm specifying which parameter for the confidence interval
#' @param level level of confidence interval, default = 0.95
#' @param parallel estimation being calculated parallely, default = FALSE
#' @param threads number of workers using to run the estimation, default = 4
#' @name confint.blblm
#' @export
#' @method confint blblm
confint.blblm <- function(object, parm = NULL, level = 0.95, parallel = FALSE, threads = 4, ...) {
  if (is.null(parm)) {
    parm <- attr(terms(object$formula), "term.labels")
  }
  alpha <- 1 - level
  est <- object$estimates
  out <- map_rbind(parm, function(p) {
    map_mean(est, ~ map_dbl(., list("coef", p)) %>% quantile(c(alpha / 2, 1 - alpha / 2)), parallel, threads)
  })
  if (is.vector(out)) {
    out <- as.matrix(t(out))
  }
  dimnames(out)[[1]] <- parm
  out
}

#' @param object blblm object
#' @param new_data new dataset
#' @param confidence confidence interval, default = FALSE
#' @param level level of prediction interval, default = 0.95
#' @param parallel estimation being calculated parallely, default = FALSE
#' @param threads number of workers using to run the estimation, default = 4
#' @name predict.blblm
#' @export
#' @method predict blblm
predict.blblm <- function(object, new_data, confidence = FALSE, level = 0.95, parallel = FALSE, threads = 4, ...) {
  est <- object$estimates
  X <- model.matrix(reformulate(attr(terms(object$formula), "term.labels")), new_data)
  if (confidence) {
    map_mean(est, ~ map_cbind(., ~ X %*% .$coef) %>%
      apply(1, mean_lwr_upr, level = level) %>%
      t(), parallel, threads)
  } else {
    map_mean(est, ~ map_cbind(., ~ X %*% .$coef) %>% rowMeans(), parallel, threads)
  }
}


mean_lwr_upr <- function(x, level = 0.95) {
  alpha <- 1 - level
  c(fit = mean(x), quantile(x, c(alpha / 2, 1 - alpha / 2)) %>% set_names(c("lwr", "upr")))
}

map_mean <- function(.x, .f, parallel, threads, ...) {
  if (parallel){
    suppressWarnings(plan(multiprocess, workers = threads))
    options(future.rng.onMisuse = "ignore")
    (future_map(.x, .f, ...) %>% reduce(`+`)) / length(.x)
  }else{
    (map(.x, .f, ...) %>% reduce(`+`)) / length(.x)
  }
}

map_cbind <- function(.x, .f, ...) {
  map(.x, .f, ...) %>% reduce(cbind)
}

map_rbind <- function(.x, .f, ...) {
  map(.x, .f, ...) %>% reduce(rbind)
}
