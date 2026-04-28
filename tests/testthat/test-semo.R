test_that("semo recovers true betas on clean bivariate data", {
  set.seed(42)
  K <- 80
  b_exp  <- rnorm(K, 0, 0.1);  se_exp <- rep(0.01, K)
  theta1 <- 0.3; theta2 <- 0.5
  b_out_1 <- theta1 * b_exp + rnorm(K, 0, 0.01)
  b_out_2 <- theta2 * b_exp + rnorm(K, 0, 0.01)
  b_out  <- cbind(b_out_1, b_out_2)
  se_out <- matrix(0.01, K, 2)
  
  fit <- fusiomr(b_exp, se_exp, b_out, se_out,
                 model = "semo",
                 control = parameter_control(niter = 3000),
                 verbose = FALSE)
  
  expect_type(fit, "list")
  expect_length(fit$est, 2)
  expect_length(fit$se, 2)
  expect_true(all(is.finite(fit$est)))
  expect_true(abs(fit$est[1] - theta1) < 5 * fit$se[1])
  expect_true(abs(fit$est[2] - theta2) < 5 * fit$se[2])
  # CI coverage
  expect_true(fit$ci[1, 1] < theta1 && fit$ci[1, 2] > theta1)
  expect_true(fit$ci[2, 1] < theta2 && fit$ci[2, 2] > theta2)
})