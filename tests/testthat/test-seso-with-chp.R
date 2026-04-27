test_that("seso_with_chp recovers true beta on clean data (no CHP)", {
  set.seed(42)
  K <- 50
  b_exp  <- rnorm(K, 0, 0.1);  se_exp <- rep(0.01, K)
  theta  <- 0.3
  b_out  <- theta * b_exp + rnorm(K, 0, 0.01)
  se_out <- rep(0.01, K)
  
  fit <- fusiomr(b_exp, se_exp, b_out, se_out,
                 model = "seso_with_chp",
                 control = parameter_control(niter = 3000),
                 verbose = FALSE)
  
  expect_type(fit, "list")
  expect_true(all(c("est", "se", "pval", "ci", "q", "n_iv") %in% names(fit)))
  expect_true(is.finite(fit$est))
  expect_true(fit$se > 0)
  expect_true(abs(fit$est - theta) < 5 * fit$se)
  expect_true(fit$ci[1] < theta && fit$ci[2] > theta)
})