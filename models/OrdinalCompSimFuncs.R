
################################################################################

# Cohen et al. 2026                                                 
# Functions to create simulated datasets

################################################################################


## load in model functions ##

library(dplyr)

source('models/OrdinalCompModelFuncs.R')

gen_sim_data = function(L=5,  # number of categories
                        N_is=500,  # number of in-sample seqs
                        N_oos=1000,  # number of out-of-sample seqs
                        zeta = 1e-12,  # compositional rounding,
                        covariates=NULL, 
                        seed=round(runif(1, min=1, max=1e9))) {
  
  # set seed for reproducibility
  set.seed(seed)
  
  # data size #
  N = N_is + N_oos
  r = sample(c(1, 3, 5, 10), N, replace=T)  # number of images in each seq
  R = sum(r)  # total number of images
  seq_ID = rep_vec(1:N, r)
  
  # covariates #
  
  if (!is.null(covariates)) {
    covariates = covariates %>% sample_n(N, replace=F)
    # note: if sample is too small, scaled binary params cause errors
  } else {
    covariates = rnorm(N*6, 0, 1) |> matrix(ncol=6)
    colnames(covariates) = paste('Var', 1:ncol(covariates), sep='')
    covariates = cbind(covariates, is_day=rbinom(N, 1, 0.5))
    covariates = as.data.frame(covariates)
  }
  
  # true data process #
  
  # generate X
  X = covariates %>% select(-is_day)
  X = scale(X)
  
  # max probability Y is in each central category
  Y_probs = c(0.2, 0.5, 0.3)
  
  # Z cutoffs theta
  theta_tilde = sapply(Y_probs, \(x) default_theta_fn(x, 3)$theta_tilde)
  theta = theta_fn(theta_tilde)
  
  # true category coefficients beta0, beta
  beta0 = median(theta[-c(1,L+1)]) + 0.2  
  # beta = sample(c(-1, -0.5, 0, 0, 0.5, 1), ncol(X), replace=F)
  beta = c(0.2, -0.3,  0.0,  0.0,  0.2, -0.3)
  names(beta) = colnames(X)
  
  
  # manual annotation process #
  
  A = 3  # number of annotators
  
  # prob Z=col given Y=row
  Z_probs = matrix(c(0.95, 0.048, 0.001, 0.0005, 0.0005,
                     0.025, 0.95, 0.024, 0.0005, 0.0005,
                     0.0005, 0.0245, 0.95, 0.0245, 0.0005,
                     0.0005, 0.0005, 0.024, 0.95, 0.025,
                     0.0005, 0.0005, 0.001, 0.048, 0.95),
                   nrow=L, byrow=T)
  
  nu_tilde = -qnorm(Z_probs[,1])  # center of annotation means (just used as a helper here)
  nu_var = 0.1  # how closely nu holds to nu_tilde
  
  # ordinal response cutoffs phi
  phi = matrix(NA, nrow=L, ncol=L+1)
  phi[,1] = -Inf
  phi[,2] = 0
  phi[,L+1] = Inf
  for (l in 3:L)
    phi[,l] = qnorm(Z_probs[,l-1]+rowSums(matrix(Z_probs[,1:(l-2)], nrow=L))) + nu_tilde
  
  phi_tilde = sapply(3:L, \(l) log(phi[,l]-phi[,l-1]))
  
  
  # AI prediction data process #
  
  # generate U
  
  U = covariates[seq_ID,] %>%
    select(is_day)
  U = scale(U)
  
  # compositional response means alpha (mean is row l given Y=l) 
  alpha = diag(5) + 0.1 - 0.5*diag(5)

  # precision coefficients omega0, omega
  omega0 = 0
  omega = 1
  names(omega) = colnames(U)
  
  # generate data #
  
  # full dataset for sampling
  data_extra = gen_data(N=N, r=r, L=L, N_oos=N_oos,
                        Z_labelled=T, C_labelled=T, 
                        beta0=beta0, beta=beta, 
                        theta_tilde=theta_tilde, phi_tilde=phi_tilde,
                        nu_tilde=nu_tilde, nu_var=nu_var, A=A, 
                        alpha=alpha, omega0=omega0, omega=omega, 
                        X=X, U=U,
                        zeta=zeta,
                        seed=seed)
  
  # remove out-of-sample labels
  data_extra$data_oos$Z = NULL
  data_extra$data_oos$C = NULL
  data_extra$data_oos$C_tilde = NULL
  
  return(data_extra)
  
}



convert_sim_data = function(data_extra_full, setting,
                            max_Z_labelled=data_extra_full$data$Z_labelled, 
                            max_C_labelled=data_extra_full$data$C_labelled, 
                            thresholds=c(0.9)) {
  
  data_full = data_extra_full$data
  extra_full = data_extra_full$extra
  
  new_Z_labelled = rep(F, data_full$R)  # ordinal images kept
  new_C_labelled = rep(F, data_full$R)  # compositional images kept
  to_max = F  # whether to convert to max value
  to_linear = F  # whether to convert to linear value
  C_threshold = 0  # theshold for AI confidences
  
  
  if (setting==1) {  # setting 1 ~ linear
    to_linear = T
    new_Z_labelled = max_Z_labelled
    new_C_labelled = max_C_labelled
    C_threshold = thresholds[setting]
  } else if (setting==2) {  # setting 2 ~ threshold
    to_max = T
    new_Z_labelled = max_Z_labelled
    new_C_labelled = max_C_labelled
    C_threshold = thresholds[setting-length(thresholds)]
  } else if (setting==3) {  # setting 3 ~ ordinal only
    new_Z_labelled = max_Z_labelled
  } else if (setting==4) {  # setting 4 ~ compositional only
    new_C_labelled = max_C_labelled
  } else if (setting==5) {  # setting 5 ~ full model
    new_Z_labelled = max_Z_labelled
    new_C_labelled = max_C_labelled
  }  
  fallible = !to_max
  
  # convert data
  data_extra_conv = convert_data(data_full, extra_full, 
                                 new_Z_labelled=new_Z_labelled,
                                 new_C_labelled=new_C_labelled,
                                 to_max=to_max, to_linear=to_linear, C_threshold=C_threshold)
  
  data_extra_conv$fallible = fallible
  
  return(data_extra_conv)
  
}

