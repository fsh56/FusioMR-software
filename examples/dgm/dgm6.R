
# function to generate two exposures two outcome GWAS summary statistics
# all IVs have UHP effect and some have CHP
# UHP effect size, CHP indicator for two outcomes are correlated

# function to generate two correlated uniform using Gaussian copula
gen.gauss.cop <- function(r, n){
rho <- 2 * sin(r * pi/6)        # Pearson correlation
P <- toeplitz(c(1, rho))        # Correlation matrix
d <- nrow(P)                    # Dimension
## Generate sample
U <- pnorm(matrix(rnorm(n*d), ncol = d) %*% chol(P))
return(U)
}

# function to map correlation to bi-Binomial prob.
bern_prob <- function(rho, p, q) {
p00 = rho * sqrt(p*q*(1-p)*(1-q)) + (1-p)*(1-q)
p10 = 1 - q - p00
p01 = 1 - p - p00
p11 = p + q + p00 -  1
return(c(p00, p01, p10, p11))
}

dgm_individual_6 = function (m, nx1, nx2, ny1, ny2, a_gamma1, b_gamma1, 
a_gamma2, b_gamma2, rho_gamma, a_f, b_f, a_alpha1, b_alpha1, a_alpha2, b_alpha2, rho_alpha, 
a_phi1, b_phi1, a_phi2, b_phi2, rho_eta, theta1, theta2, q_uhp1=1, q_uhp2=1, 
q_chp1=0.1, q_chp2=0.1, beta_XU=1, beta_YU=1) {
# IV-to-exposure effect (not through U)
gamma = gen.gauss.cop(rho_gamma, m)
gamma1 = a_gamma1 + (b_gamma1 - a_gamma1) * gamma[,1]
gamma2 = a_gamma2 + (b_gamma2 - a_gamma2) * gamma[,2]

# MAF
f1x = runif(m, a_f, b_f)
f1y = runif(m, a_f, b_f)
f2x = runif(m, a_f, b_f)
f2y = runif(m, a_f, b_f)
# Genotypes in [0, 1, 2], m x n
gx1 = replicate(nx1, rbinom(n = m, size = 2, prob = f1x))
gx2 = replicate(nx2, rbinom(n = m, size = 2, prob = f2x))
gy1 = replicate(ny1, rbinom(n = m, size = 2, prob = f1y))
gy2 = replicate(ny2, rbinom(n = m, size = 2, prob = f2y))
# for each row (snp), standardize genotype
gx1 = t(apply(gx1, 1, base::scale))
gx2 = t(apply(gx2, 1, base::scale))
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

# CHP
# eta, correlated Bernoulli
bernp = bern_prob(rho_eta, q_chp1, q_chp2)
d = rmultinom(m, 1, bernp)
# print(rowSums(d)/K)
eta = t(matrix(c(0,0,0,1,1,0,1,1), nrow = 2) %*% d)
# cor(eta)
# index of invalid IVs with CHP
ind1 = which(eta[,1] == 1)
ind2 = which(eta[,2] == 1)
# CHP effect
phi1 = phi2 = rep(0, m)
phi1[ind1] = runif(length(ind1), a_phi1, b_phi1)
phi2[ind2] = runif(length(ind2), a_phi2, b_phi2)

# exposure and outcome
Ux1 = c(phi1 %*% gx1) + rnorm(nx1)
Ux2 = c(phi2 %*% gx2) + rnorm(nx2)
Uy1 = c(phi1 %*% gy1) + rnorm(ny1)
Uy2 = c(phi2 %*% gy2) + rnorm(ny2)
Xx1 = c(gamma1 %*% gx1) + beta_XU * Ux1 + rnorm(nx1)
Xx2 = c(gamma2 %*% gx2) + beta_XU * Ux2 + rnorm(nx2)
Xy1 = c(gamma1 %*% gy1) + beta_XU * Uy1 + rnorm(ny1)
Xy2 = c(gamma2 %*% gy2) + beta_XU * Uy2 + rnorm(ny2)
Y1 = c(alpha1 %*% gy1) + theta1 * Xy1 + beta_YU * Uy1 + rnorm(ny1)
Y2 = c(alpha2 %*% gy2) + theta2 * Xy2 + beta_YU * Uy2 + rnorm(ny2)
# variance component 
h2_exposure1 = var(c(gamma1 %*% gy1))/var(Xy1)
h2_exposure2 = var(c(gamma2 %*% gy2))/var(Xy2)
h2_causal1 = var(c(gamma1 %*% gy1 * theta1))/var(Y1)
h2_causal2 = var(c(gamma2 %*% gy2 * theta2))/var(Y2)
h2_uhp1 = var(c(alpha1 %*% gy1))/var(Y1)
h2_uhp2 = var(c(alpha2 %*% gy2))/var(Y2)
h2_chp1 = var(c(phi1 %*% gy1 * beta_YU))/var(Y1)
h2_chp2 = var(c(phi2 %*% gy2 * beta_YU))/var(Y2)

# output
return(list(X1 = Xx1, X2 = Xx2, Y1 = Y1, Y2 = Y2, gx1 = gx1, gx2 = gx2, 
gy1 = gy1, gy2 = gy2, h2_exposure1 = h2_exposure1, h2_exposure2 = h2_exposure2, 
h2_causal1 = h2_causal1, h2_causal2 = h2_causal2, h2_uhp1 = h2_uhp1, h2_uhp2 = h2_uhp2, 
h2_chp1 = h2_chp1, h2_chp2 = h2_chp2))
}

# generate summary statistics
dgm_summary_6 = function(GWAS_individual) {
gammaSS1 = fastSigLm(GWAS_individual$X1, t(GWAS_individual$gx1))
gammaSS2 = fastSigLm(GWAS_individual$X2, t(GWAS_individual$gx2))
GammaSS1 = fastSigLm(GWAS_individual$Y1, t(GWAS_individual$gy1))
GammaSS2 = fastSigLm(GWAS_individual$Y2, t(GWAS_individual$gy2))
b_exp_1 = gammaSS1$coef; se_exp_1 = gammaSS1$std
b_exp_2 = gammaSS2$coef; se_exp_2 = gammaSS2$std
b_out_1 = GammaSS1$coef; se_out_1 = GammaSS1$std  
b_out_2 = GammaSS2$coef; se_out_2 = GammaSS2$std  
return(list(
b_exp_1 = b_exp_1, se_exp_1 = se_exp_1, b_exp_2 = b_exp_2, se_exp_2 = se_exp_2, 
b_out_1 = b_out_1, se_out_1 = se_out_1, b_out_2 = b_out_2, se_out_2 = se_out_2,
h2_exposure1 = GWAS_individual$h2_exposure1, h2_exposure2 = GWAS_individual$h2_exposure2,
h2_causal1 = GWAS_individual$h2_causal1, h2_causal2 = GWAS_individual$h2_causal2,
h2_uhp1 = GWAS_individual$h2_uhp1, h2_uhp2 = GWAS_individual$h2_uhp2,
h2_chp1 = GWAS_individual$h2_chp1, h2_chp2 = GWAS_individual$h2_chp2))
}

# wrap up
dgm6 = function(m, nx1, nx2, ny1, ny2, a_gamma1, b_gamma1, 
a_gamma2, b_gamma2, rho_gamma, a_f, b_f, a_alpha1, b_alpha1, a_alpha2, b_alpha2, rho_alpha, 
a_phi1, b_phi1, a_phi2, b_phi2, rho_eta, theta1, theta2, q_uhp1=1, q_uhp2=1, 
q_chp1=0.1, q_chp2=0.1, beta_XU=1, beta_YU=1) {
GWAS_individual = dgm_individual_6(m, nx1, nx2, ny1, ny2, a_gamma1, b_gamma1, 
a_gamma2, b_gamma2, rho_gamma, a_f, b_f, a_alpha1, b_alpha1, a_alpha2, b_alpha2, rho_alpha, 
a_phi1, b_phi1, a_phi2, b_phi2, rho_eta, theta1, theta2, q_uhp1, q_uhp2, 
q_chp1, q_chp2, beta_XU, beta_YU) 
GWAS_summary = dgm_summary_6(GWAS_individual)
return(GWAS_summary)
}


