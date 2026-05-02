
# function to generate a single population GWAS summary statistics
# all generated IVs are used. IV strength is controlled by a two-disjoint-interval uniform
# allow proportional CHP
# allow correlated IV. LD correlation = AR(1)

ar1_cor <- function(n, rho) {
exponent = abs(matrix(1:n - 1, nrow = n, ncol = n, byrow = TRUE) - (1:n - 1))
rho^exponent
}

# Gaussian-copula AR(1) genotypes with a provided LD matrix (Sigma)
simulate_genotypes_from_ld <- function(m, n, maf, Sigma) {
Z1 = mvrnorm(n = n, mu = rep(0, m), Sigma = Sigma)  # n x m
Z2 = mvrnorm(n = n, mu = rep(0, m), Sigma = Sigma)
thr = qnorm(maf)
H1 = sweep(Z1, 2, thr, "<") * 1L
H2 = sweep(Z2, 2, thr, "<") * 1L
t(H1 + H2)  # m x n (0/1/2)
}


dgm_individual_4_ld = function (m, nx, ny, a_gamma, b_gamma,
a_f, b_f, a_alpha, b_alpha, a_phi, b_phi, theta, q_uhp, q_chp, beta_XU = 1, beta_YU = 1, rho_ld = 0) {
# IV-to-exposure effect (not through U)
gamma = runif(m, a_gamma, b_gamma)
# MAF
f = runif(m, a_f, b_f)
# correlated genotypes, in [0, 1, 2], m x n, LD matrix (shared by gx & gy)
if (abs(rho_ld) < 1e-12) {
R_ld = diag(m)
gx = replicate(nx, rbinom(n = m, size = 2, prob = f))
gy = replicate(ny, rbinom(n = m, size = 2, prob = f))
} else {
R_ld = ar1_cor(m, rho_ld)
gx = simulate_genotypes_from_ld(m, nx, maf = f, Sigma = R_ld)
gy = simulate_genotypes_from_ld(m, ny, maf = f, Sigma = R_ld)
}

# for each row (snp), standardize genotype
gx = t(apply(gx, 1, base::scale))
gy = t(apply(gy, 1, base::scale))
# a proportion of snps have UHP effect (q_uhp)
alpha = runif(m, a_alpha, b_alpha)
valid_ids = sample(1:m, floor((1-q_uhp)*m))
alpha[valid_ids] = 0
# a proportion of snps have CHP effect (q_chp)
eta = rbinom(m, 1, q_chp) 
ind_chp = which(eta == 1) 
phi = rep(0, m)
phi[ind_chp] = runif(length(ind_chp), a_phi, b_phi)

# exposure and outcome
Ux = c(phi %*% gx) + rnorm(nx)
Uy = c(phi %*% gy) + rnorm(ny)
Xx = c(gamma %*% gx) + beta_XU * Ux + rnorm(nx)
Xy = c(gamma %*% gy) + beta_XU * Uy + rnorm(ny)
Y  = c(alpha %*% gy) + theta * Xy + beta_YU * Uy + rnorm(ny)

h2_exposure = var(c(gamma %*% gy))/var(Xy)
h2_causal = var(c(gamma %*% gy * theta))/var(Y)
h2_uhp = var(c(alpha %*% gy))/var(Y)
h2_chp = var(c(phi %*% gy * beta_YU))/var(Y)
# output
return(list(X = Xx, Y = Y, gx = gx, gy = gy, h2_exposure = h2_exposure, 
h2_causal = h2_causal, h2_uhp = h2_uhp, h2_chp = h2_chp, R_ld = R_ld))
}

# generate summary statistics
dgm_summary_4 = function(GWAS_individual) {
gammaSS = fastSigLm(GWAS_individual$X, t(GWAS_individual$gx))
GammaSS = fastSigLm(GWAS_individual$Y, t(GWAS_individual$gy))
b_exp = gammaSS$coef; se_exp = gammaSS$std
b_out = GammaSS$coef; se_out = GammaSS$std  
return(list(b_exp = b_exp, se_exp = se_exp, b_out = b_out, se_out = se_out,
h2_exposure = GWAS_individual$h2_exposure, h2_causal = GWAS_individual$h2_causal,
h2_uhp = GWAS_individual$h2_uhp, h2_chp = GWAS_individual$h2_chp))
}

# wrap up
dgm4_ld = function(m, nx, ny, a_gamma, b_gamma, 
a_f, b_f, a_alpha, b_alpha, a_phi, b_phi, theta, q_uhp, q_chp, beta_XU = 1, beta_YU = 1, rho_ld = 0) {
GWAS_individual = dgm_individual_4_ld(m, nx, ny, a_gamma, b_gamma, a_f, b_f, a_alpha, b_alpha, a_phi, b_phi, theta, q_uhp, q_chp, beta_XU = 1, beta_YU = 1, rho_ld) 
GWAS_summary = dgm_summary_4(GWAS_individual)
GWAS_summary$R_ld = GWAS_individual$R_ld
#cor(t(GWAS_individual$gy[1:5,])); GWAS_individual$R_ld[1:5,1:5]; cor(t(GWAS_individual$gx[1:5,]))
return(GWAS_summary)
}

