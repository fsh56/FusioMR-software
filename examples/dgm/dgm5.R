
# function to generate single exposure, two outcome GWAS summary statistics
# all generated IVs are used. 
# UHP only
# UHP effect for two outcomes are correlated

# function to generate two correlated uniform using Gaussian copula
gen.gauss.cop <- function(r, n){
rho <- 2 * sin(r * pi/6)        # Pearson correlation
P <- toeplitz(c(1, rho))        # Correlation matrix
d <- nrow(P)                    # Dimension
## Generate sample
U <- pnorm(matrix(rnorm(n*d), ncol = d) %*% chol(P))
return(U)
}

dgm_individual_5 = function (m, nx, ny1, ny2, a_gamma, b_gamma,
a_f, b_f, a_alpha1, b_alpha1, a_alpha2, b_alpha2, rho_alpha, theta1, theta2, 
q_uhp1 = 1, q_uhp2 = 1) {
# IV-to-exposure effect (not through U)
gamma = runif(m, a_gamma, b_gamma)

# MAF
f = runif(m, a_f, b_f)
# Genotypes in [0, 1, 2], m x n
gx = replicate(nx, rbinom(n = m, size = 2, prob = f))
gy1 = replicate(ny1, rbinom(n = m, size = 2, prob = f))
gy2 = replicate(ny2, rbinom(n = m, size = 2, prob = f))
# for each row (snp), standardize genotype
gx = t(apply(gx, 1, base::scale))
gy1 = t(apply(gy1, 1, base::scale))
gy2 = t(apply(gy2, 1, base::scale))

# a proportion of snps have UHP effect (q_uhp)
alpha = gen.gauss.cop(rho_alpha, m)
alpha[,1] = a_alpha1 + (b_alpha1 - a_alpha1) * alpha[,1]
alpha[,2] = a_alpha2 + (b_alpha2 - a_alpha2) * alpha[,2]
valid_ids1 = sample(1:m, floor((1-q_uhp1)*m))
valid_ids2 = sample(1:m, floor((1-q_uhp2)*m))
alpha[valid_ids1, 1] = 0
alpha[valid_ids2, 2] = 0
alpha1 = alpha[,1]
alpha2 = alpha[,2]

# exposure and outcome
phi = rep(0, m)
Ux = c(phi %*% gx) + rnorm(nx)
Uy1 = c(phi %*% gy1) + rnorm(ny1)
Uy2 = c(phi %*% gy2) + rnorm(ny2)
beta_XU = beta_YU = 1
Xx = c(gamma %*% gx) + beta_XU * Ux + rnorm(nx)
Xy1 = c(gamma %*% gy1) + beta_XU * Uy1 + rnorm(ny1)
Xy2 = c(gamma %*% gy2) + beta_XU * Uy2 + rnorm(ny2)
Y1 = c(alpha1 %*% gy1) + theta1 * Xy1 + beta_YU * Uy1 + rnorm(ny1)
Y2 = c(alpha2 %*% gy2) + theta2 * Xy2 + beta_YU * Uy2 + rnorm(ny2)

h2_exposure1 = var(c(gamma %*% gy1))/var(Xy1)
h2_exposure2 = var(c(gamma %*% gy2))/var(Xy2)
h2_causal1 = var(c(gamma %*% gy1 * theta1))/var(Y1)
h2_causal2 = var(c(gamma %*% gy2 * theta2))/var(Y2)
h2_uhp1 = var(c(alpha1 %*% gy1))/var(Y1)
h2_uhp2 = var(c(alpha2 %*% gy2))/var(Y2)

# output
return(list(X = Xx, Y1 = Y1, Y2 = Y2, gx = gx, gy1 = gy1, gy2 = gy2, 
h2_exposure1 = h2_exposure1, h2_exposure2 = h2_exposure2, 
h2_causal1 = h2_causal1, h2_causal2 = h2_causal2, 
h2_uhp1 = h2_uhp1, h2_uhp2 = h2_uhp2))
}

# generate summary statistics
dgm_summary_5 = function(GWAS_individual) {
gammaSS = fastSigLm(GWAS_individual$X, t(GWAS_individual$gx))
GammaSS1 = fastSigLm(GWAS_individual$Y1, t(GWAS_individual$gy1))
GammaSS2 = fastSigLm(GWAS_individual$Y2, t(GWAS_individual$gy2))
b_exp = gammaSS$coef; se_exp = gammaSS$std
b_out_1 = GammaSS1$coef; se_out_1 = GammaSS1$std  
b_out_2 = GammaSS2$coef; se_out_2 = GammaSS2$std  
return(list(b_exp = b_exp, se_exp = se_exp, b_out_1 = b_out_1, se_out_1 = se_out_1,
b_out_2 = b_out_2, se_out_2 = se_out_2,
h2_exposure1 = GWAS_individual$h2_exposure1, h2_exposure2 = GWAS_individual$h2_exposure2,
h2_causal1 = GWAS_individual$h2_causal1, h2_causal2 = GWAS_individual$h2_causal2,
h2_uhp1 = GWAS_individual$h2_uhp1, h2_uhp2 = GWAS_individual$h2_uhp2))
}

# wrap up
dgm5 = function(m, nx, ny1, ny2, a_gamma, b_gamma,
a_f, b_f, a_alpha1, b_alpha1, a_alpha2, b_alpha2, rho_theta, theta1, theta2, 
q_uhp1 = 1, q_uhp2 = 1) {
GWAS_individual = dgm_individual_5(m, nx, ny1, ny2, a_gamma, b_gamma,
a_f, b_f, a_alpha1, b_alpha1, a_alpha2, b_alpha2, rho_theta, theta1, theta2, 
q_uhp1 = 1, q_uhp2 = 1) 
GWAS_summary = dgm_summary_5(GWAS_individual)
return(GWAS_summary)
}


