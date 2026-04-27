test_that("seso_uhp_only recovers true beta on clean data", {
  set.seed(42)
  K <- 80
  b_exp  <- rnorm(K, 0, 0.1)
  se_exp <- rep(0.01, K)
  theta  <- 0.3
  b_out  <- theta * b_exp + rnorm(K, 0, 0.01)
  se_out <- rep(0.01, K)

  fit <- fusiomr(b_exp, se_exp, b_out, se_out,
                 model = "seso_uhp_only",
                 verbose = FALSE)

  expect_type(fit, "list")
  expect_true(all(c("est", "se", "pval", "ci", "n_iv") %in% names(fit)))
  expect_true(is.finite(fit$est))
  expect_true(fit$se > 0)
  # estimate within a few SEs of the truth
  expect_true(abs(fit$est - theta) < 5 * fit$se)
  # 95% CI covers truth
  expect_true(fit$ci[1] < theta && fit$ci[2] > theta)
})

test_that("fusiomr input validation works", {
  expect_error(fusiomr("a", 1, 1, 1, model = "seso_uhp_only"),
               "numeric")
  expect_error(fusiomr(c(1, 2), c(1), c(1, 2), c(1, 2),
                       model = "seso_uhp_only"),
               "same length")
  expect_error(fusiomr(c(1, 2), c(1, 2), c(1, 2), c(-1, 2),
                       model = "seso_uhp_only", control = parameter_control(niter = 500)),
               "positive")
})

test_that("parameter_control returns a well-formed list", {
  ctrl <- parameter_control()
  expect_type(ctrl, "list")
  expect_true(all(c("niter", "burnin_prop",
                    "c_gamma", "c_theta", "kappa_gamma", "kappa_theta",
                    "rho_ov", "z_thresh") %in% names(ctrl)))
  expect_equal(ctrl$c_gamma, 0.5)
  expect_equal(ctrl$niter, 20000)
})
