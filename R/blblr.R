#' @import purrr
#' @import furrr
#' @import stats
#' @import tidyverse
#' @import utils
#' @importFrom magrittr %>%
#' @details
#' Linear Regression with Little Bag of Bootstraps
"_PACKAGE"

## quiets concerns of R CMD check re: the .'s that appear in pipelines
# from https://github.com/jennybc/googlesheets/blob/master/R/googlesheets.R
utils::globalVariables(c("."))

#' @param formula logistic regression formula
#' @param data dataset
#' @param m number of split, default = 10
#' @param B number of bootstrap, default = 5000
#' @param parallel estimation being calculated parallely, default = FALSE
#' @param threads number of workers using to run the estimation, default = 4
#' @export
#' Logistic regression version of bag of little bootstraps
blblr <- function(formula, data, m = 10, B = 5000, parallel = FALSE, threads = 4){
  if (parallel){
    suppressWarnings(future::plan(multiprocess, workers = threads))
    data_list <- split_data(data, m)
    estimates <- future_map(data_list,~lr_each_subsample(formula = formula, data = ., n = nrow(data), B = B))
  }else{
    data_list <- split_data(data, m)
    estimates <- map(data_list, ~lr_each_subsample(formula = formula, data = ., n = nrow(data), B = B))
  }
  res <- list(estimates = estimates, formula = formula)
  class(res) <- "blblr"
  invisible(res)
}

#' @param data dataset
#' @param m number of splits
#' split data into m parts of approximated equal sizes
split_data <- function(data, m) {
  idx <- sample.int(m, nrow(data), replace = TRUE)
  data %>% split(idx)
}

#' compute the estimates
lr_each_subsample <- function(formula, data, n, B) {
  freqs <- rmultinom(1, n, rep(1, nrow(data)))
  replicate(B, lr1(formula, data, freqs), simplify = FALSE)
}

#' estimate the regression estimates based on given the number of repetitions
lr1 <- function(formula, data, freqs){
  environment(formula) <- environment()
  fit <- glm(formula, family="binomial", data, weights = freqs, maxit=100)
  list(coef = blbcoef(fit))
}

#' @param X a design matrix of dimension n * p
#' @param y a vector of observations of length n, or a matrix with n rows
#' @param n number of rows in each subset
#' compute the regression estimates for a blb dataset
lr1 <- function(X, y, n){
  freqs <- as.vector(rmultinom(1, n, rep(1, nrow(X))))
  fit <- glm.fit(X, y, freqs)
  list(coef = blbcoef(fit))
}

#' @param x dataset for blblr model
#' @export
#' @method print blblr
print.blblr <- function(x, ...) {
  cat("blblrodel:", capture.output(x$formula))
  cat("\n")
}

#' @param object blblr model
#' @param parallel calculate the mean parallely, default = FALSE
#' @param threads number of workers using to run the calculation, default = 4
#' @export
#' @method coef blblr
coef.blblr <- function(object, parallel = FALSE, threads = 4,...) {
  est <- object$estimates
  map_mean(est, ~ map_cbind(., "coef") %>% rowMeans(), parallel, threads)
}


#' @param object blblr object
#' @param parm specifying which parameter for the confidence interval
#' @param level level of confidence interval, default = 0.95
#' @param parallel calculate the mean parallelly, default = FALSE
#' @param threads number of workers using to run the calculation, default = 4
#' @name confint.blblr
#' @export
#' @method confint blblr
confint.blblr <- function(object, parm = NULL, level = 0.95, parallel = FALSE, threads = 4, ...) {
  if (is.null(parm)) {
    parm <- attr(terms(object$formula), "term.labels")
  }
  alpha <- 1 - level
  est <- object$estimates
  out <- map_rbind(parm, function(p) {
    map_mean(est, ~ map_dbl(., list("coef", p)) %>% quantile(c(alpha / 2, 1 - alpha / 2))
             , parallel, threads)
  })
  if (is.vector(out)) {
    out <- as.matrix(t(out))
  }
  dimnames(out)[[1]] <- parm
  out
}

#' @param object blblr object
#' @param new_data new dataset
#' @param confidence confidence interval, default = FALSE; TRUE --> prediction interval
#' @param level level of prediction interval, default = 0.95
#' @param parallel calculate the mean parallelly, default = FALSE
#' @param threads number of workers using to run the calculation, default = 4
#' @name predict.blblr
#' @export
#' @method predict blblr
predict.blblr <- function(object, new_data, confidence = FALSE, level = 0.95, parallel = FALSE, threads = 4, ...) {
  est <- object$estimates
  X <- model.matrix(reformulate(attr(terms(object$formula), "term.labels")), new_data)
  if (confidence) {
    map_mean(est, ~ map_cbind(., ~ X %*% .$coef) %>%
               apply(1, mean_lwr_upr, level = level) %>%
               t(), parallel, threads)
  } else {
    map_mean(est, ~ map_cbind(., ~ X %*% .$coef) %>% rowMeans())
  }
}


mean_lwr_upr <- function(x, level = 0.95) {
  alpha <- 1 - level
  c(fit = mean(x), quantile(x, c(alpha / 2, 1 - alpha / 2)) %>% set_names(c("lwr", "upr")))
}

map_mean <- function(.x, .f, parallel, threads, ...) {
  if (parallel){
    suppressWarnings(future::plan(multiprocess, workers = threads))
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