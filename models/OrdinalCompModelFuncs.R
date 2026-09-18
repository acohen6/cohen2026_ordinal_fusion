
################################################################################

# Cohen et al. 2026                                                 
# Functions to support model fitting

################################################################################


## imports ##

library(rjags, coda)  # for linear fit

library(LaplacesDemon)  # rdirichlet


## default setting functions ##

# parameters that can be sampled given available data
default_sampling_params_fn = function(data,
                                      sampling_params=NA,
                                      fixed_params=NA,
                                      fallible=T) {
  # default params if none specified
  default_sampling_params = c('Y', 'beta0', 'beta', 'theta_tilde', 
                              'nu', 'nu_tilde', 'phi_tilde', 
                              'alpha', 'omega0', 'omega')
  if (sum(!is.na(data$Z))==0)
    fixed_params = union(fixed_params, c('nu', 'nu_tilde', 'phi_tilde'))
  if (sum(!is.na(data$C))==0)
    fixed_params = union(fixed_params, c('alpha', 'omega0', 'omega'))
  if (!fallible)
    fixed_params = union(fixed_params, c('nu', 'nu_tilde', 'phi_tilde'))
  if (data$p==0)
    fixed_params = union(fixed_params, c('beta'))
  if (data$q==0)
    fixed_params = union(fixed_params, c('omega'))
  
  # check only sampling valid params
  if (any(is.na(sampling_params))) {
    sampling_params = default_sampling_params
  } else {
    sampling_params = intersect(sampling_params, default_sampling_params)
  }
  
  sampling_params = sampling_params[!sampling_params %in% fixed_params]
  
  return(sampling_params)
}

# names for parameters given size of data
default_param_names_fn = function(data_size) {
  
  # names for sampling parameters (match size of parameters sampled)
  sampling_param_names = list(Y=paste('Y', 1:data_size$N, sep='_'),
                              beta0='beta0',
                              beta=paste('beta', 1:data_size$p, sep='_'),
                              theta_tilde=paste('theta_tilde', 2:(data_size$L-1), sep='_'),
                              nu=paste(rep(paste('nu', 1:data_size$A, sep='_'), data_size$L),
                                       rep(1:data_size$L, each=data_size$A), sep='_'),
                              nu_tilde=paste('nu_tilde', 1:data_size$L, sep='_'),
                              phi_tilde=paste(rep(paste('phi_tilde', 1:data_size$L, sep='_'), data_size$L-2),
                                             rep(2:(data_size$L-1), each=data_size$L-2), sep='_'),  
                              alpha=paste(rep(paste('alpha', 1:data_size$L, sep='_'), data_size$L), 
                                          rep(1:data_size$L, each=data_size$L), sep='_'),
                              omega0='omega0',
                              omega=paste('omega', 1:data_size$q, sep='_'))
  
  # names for acceptance rates (match size of proposal blocks)
  acceptance_param_names = list(beta0='beta0',
                                beta='beta',
                                theta_tilde='theta_tilde',
                                nu=paste(rep(paste('nu', 1:data_size$A, sep='_'), data_size$L),
                                         rep(1:data_size$L, each=data_size$A), sep='_'),
                                nu_tilde=paste('nu_tilde', 1:data_size$L, sep='_'),
                                phi_tilde=paste('phi_tilde', 1:data_size$L, sep='_'),
                                alpha=paste('alpha', 1:data_size$L, sep='_'),
                                omega0='omega0',
                                omega='omega')
  
  return(list(sampling_param_names=sampling_param_names, 
              acceptance_param_names=acceptance_param_names))
}

# setup for equally spaced theta with given max prob between cutoffs
default_theta_fn = function(max_Y_prob, L) {
  theta0 = 2*qnorm((max_Y_prob+1)/2)  # space between cutoffs
  theta_tilde = rep(log(theta0), L-2)
  theta = theta_fn(theta_tilde)
  beta0 = theta0/2*(L-2)
  return(list(theta0=theta0, theta_tilde=theta_tilde, theta=theta, beta0=beta0))
}

# setup for equally spaced phi with given max prob between cutoffs
default_phi_fn = function(max_Z_prob, L, A=1) {
  Z_probs = (max_Z_prob-(1-max_Z_prob)/(L-1))*diag(L)+(1-max_Z_prob)/(L-1)
  nu_tilde = -qnorm(Z_probs[,1])  
  nu = matrix(nu_tilde, nrow=A, ncol=L, byrow=T)
  
  phi = matrix(NA, nrow=L, ncol=L+1)
  phi[,1] = -Inf
  phi[,2] = 0
  phi[,L+1] = Inf
  for (l in 3:L)
    phi[,l] = qnorm(Z_probs[,l-1]+rowSums(matrix(Z_probs[,1:(l-2)], nrow=L))) + nu_tilde
  
  phi_tilde = sapply(3:L, \(l) log(phi[,l]-phi[,l-1]))
  
  return(list(nu=nu, nu_tilde=nu_tilde, phi_tilde=phi_tilde, phi=phi))
}

# default initial parameters
default_params0_fn = function(data_size, 
                              params0=NA,  # optional named list 
                              default=c()  # parameters to override info in params0
                              ) {
  
  # default initial parameters
  default_Y = rcat(data_size$N, rep(1/data_size$L, data_size$L))
  default_theta = default_theta_fn(0.5, data_size$L)
  default_phi = default_phi_fn(0.95, data_size$L, data_size$A)
  
  default_params0 = list(Y=default_Y,
                         beta0=default_theta$beta0,
                         beta=rep(0, data_size$p),
                         theta_tilde=default_theta$theta_tilde,
                         nu_tilde=default_phi$nu_tilde,
                         nu=default_phi$nu,
                         phi_tilde=default_phi$phi_tilde,
                         alpha=0.7*diag(data_size$L)+0.3/(data_size$L-1)*(1-diag(data_size$L)),
                         omega0=1,
                         omega=rep(0, data_size$q))

  # fill in params not mentioned in params0 with default settings  
  if (!'list' %in% class(params0))
    params0 = list()
  for (.param in union(names(params0), names(default_params0))) {
    if (! .param %in% names(default_params0))
      params0[[.param]] = NULL
    if (all(is.null(params0[[.param]])) | .param %in% default) {
      params0[[.param]] = default_params0[[.param]]
    }
  }
  
  # ensure nu_tilde is valid
  nu_tilde_bds = nu_tilde_bds_fn(phi_fn(params0$phi))
  for (l in 1:data_size$L) {
    if (nu_tilde_bds[l,1] > params0$nu_tilde[l] | params0$nu_tilde[l] > nu_tilde_bds[l,2])
      params0$nu_tilde[l] = mean(nu_tilde_bds[l,])
  }
  
  # return initial parameters
  return(params0)
}

# default hyperparameters (prior and initial proposal)
default_hp_fn = function(data, 
                         prior0=NA, prop0=NA,
                         sampling_params=c('Y', 'beta0', 'beta', 'theta_tilde', 
                                           'nu', 'nu_tilde', 'phi_tilde', 
                                           'alpha', 'omega0', 'omega')) {
  
  hp = list()
  
  # prior hyperparameters #
  
  # defaults
  hp$prior = prior0
  prior_default = list(max_Y_prob=1/data$L,
                       beta0_var=10,
                       beta_var=1,
                       theta_var=10,
                       max_Z_prob=0.95,
                       nu_var=0.1, 
                       phi_var=10,
                       mu_alpha=0.4,
                       alpha_precision=1,
                       omega0_mn=0,
                       omega0_var=10,
                       omega_var=10)
  
  # fill in params not mentioned in prior0 with default settings 
  if (!'list' %in% class(prior0))
    hp$prior = list()
  for (.param in names(prior_default)) {
    if (all(is.null(hp$prior[[.param]]))) {
      hp$prior[[.param]] = prior_default[[.param]]
    }
  }
  
  # proposal hyperparameters #
  
  # defaults
  hp$prop = prop0
  prop_default = list(beta0=0.2, beta=0.1, theta_tilde=0.1, 
                      nu=rep(0.1, data$A*data$L),
                      nu_tilde=rep(1, data$L),
                      phi_tilde=rep(0.5, data$L),
                      alpha=rep(1000, data$L), omega0=0.05, omega=0.05)
  
  # fill in params not mentioned in prop0 with default settings 
  if (!'list' %in% class(prop0))
    hp$prop = list()
  for (.param in sampling_params) {
    if (all(is.null(hp$prop[[.param]]))) {
      num_sd = 1
      hp$prop[[.param]] = rep(prop_default[[.param]], num_sd)
    }
  }
  
  # return hyperparameters
  return(hp)
}



## data generation functions ##

# HELPER FUNCTIONS #

# multiply matrix by a vector and add an intercept 
mu_fn = function(X, beta0, beta) {
  if (!'matrix' %in% class(X))
    X = matrix(X, ncol=length(beta))
  mu = beta0
  if (length(beta)>0)
    mu = mu + X %*% beta
  return(mu)
}

# recreate theta (vector of cutoffs for Y) from theta_tilde
theta_fn = function(theta_tilde) {
  theta = c(-Inf, 0, rep(NA, length(theta_tilde)), Inf)
  for (l in 1:length(theta_tilde)) {
    theta[l+2] = theta[l+1] + exp(theta_tilde[l])
  }
  return(theta)
}

# recreate phi (vector of cutoffs for ordinal response) from phi_tilde
phi_fn = function(phi_tilde) {
  L = nrow(phi_tilde)
  phi = matrix(c(-Inf, 0, rep(NA, L-2), Inf), nrow=L, ncol=L+1, byrow=T)
  for (l in 1:(L-2)) {
    phi[,l+2] = phi[,l+1] + exp(phi_tilde[,l])
  }
  return(phi)
}

# get bounds for nu_tilde given phi
nu_tilde_bds_fn = function(phi) {
  L = nrow(phi)
  
  nu_tilde_bds = matrix(NA, nrow=L, ncol=2)
  
  for (l in 1:L) {
    nu_tilde_bds[l,] = phi[l, l:(l+1)]
    
    if (l==1) {
      nu_tilde_bds[l,1] = -3
    } else if (l==L) {
      nu_tilde_bds[l,2] = phi[l,l] + 3
    } 
  }
  
  return(nu_tilde_bds)
}

# create precisions and concentrations for compositional response
concentration_fn = function(data, params, eval=rep(T, length(data$seq_ID))) {
  C_mns = params$alpha[params$Y[data$seq_ID][eval],]
  precision = exp(mu_fn(data$U[eval,], params$omega0, params$omega)) |> c()
  
  concentration = C_mns * precision  # multiplies precision by each columns

  return(list(precision=precision, concentration=concentration))
}

# apply adjustment to fix rounding zeros in compositional data
fix_zeros = function(C, zeta=1e-3) {
  # rows with values below cutoff
  zeros = apply(C, 1, function(c) any(c < zeta))
  zeros[is.na(zeros)] = FALSE
  
  # return if nothing to fix
  if (!any(zeros))
    return(C)
  
  # add cutoff to these rows and normalize
  C[zeros,] = C[zeros,] + zeta
  if (sum(zeros) > 1) {
    C[zeros,] = C[zeros,] / rowSums(C[zeros,])
  } else {
    C[zeros,] = C[zeros,] / sum(C[zeros,])
  }
  
  return(C)
}

# rep() with vector each for vectors
rep_vec = function(x, each) {
  sapply(seq_along(x), function(i) rep(x[i], each[i])) |> 
    unlist() |> c()
}

# colSum rows of a vector or matrix C according to the replicates vector r (ASSUMES SORTED)
r_sum = function(C, r, na.rm=F) {
  N = length(r)
  indiv_sum = NA
  
  if ('matrix' %in% class(C)) {  # for a matrix
    indiv_sum = matrix(NA, nrow=n, ncol=ncol(C))
    
    start = 1
    end = 1
    for (i in 1:n) {
      if (r[i] > 0) {
        end = start + r[i] - 1
        indiv_sum[i,] = colSums(C[start:end,], na.rm=na.rm)
        start = end + 1
      }
    }
  } else {  # for a vector
    indiv_sum = rep(0, N)
    
    start = 1
    end = 1
    for (i in 1:N) {
      if (r[i] > 0) {
        end = start + r[i] - 1
        indiv_sum[i] = sum(C[start:end], na.rm=na.rm)
        start = end + 1
      } 
    }
  }
  
  return(indiv_sum)
}


# MAIN FUNCTION #

gen_data = function(N=100, r=rep(1, N),  # size of data
                    N_oos=0,  # number of samples to separate for testing
                    L=3,  # number of categories
                    A=1,  # number of annotators
                    Z_labelled=T, C_labelled=T,  # which responses are included
                    seed=round(runif(1, min=1, max=1e9)),  # random seed 
                    
                    # detault params
                    max_Y_prob=0.5, max_Z_prob=0.9, 
                    beta0=default_theta_fn(max_Y_prob, L)$beta0,
                    beta=c(), p=length(beta),
                    theta_tilde=default_theta_fn(max_Y_prob, L)$theta_tilde, 
                    phi_tilde=default_phi_fn(max_Z_prob, L, A)$phi_tilde,
                    nu_tilde=default_phi_fn(max_Z_prob, L, A)$nu_tilde, nu_var=0.1,
                    #nu=default_phi_fn(max_Z_prob, L, A)$nu, 
                    alpha=0.9*diag(L)+0.1/(L-1)*(1-diag(L)),
                    omega0=0, omega=c(), q=length(omega),
                    X=NA, U=NA,
                    zeta=1e-8
) {
  
  # set seed for reproducibility
  set.seed(seed)
  
  ## setup ##
  
  # size of data
  if (length(r)==1)
    r = rep(r, N)
  
  R = sum(r)
  
  if (length(Z_labelled)==1)
    Z_labelled = rep(Z_labelled, R)
  if (length(C_labelled)==1)
    C_labelled = rep(C_labelled, R)
  
  seq_ID = rep_vec(1:N, r)  # sequence that each image belongs to
  
  
  # X: individual-level covariates
  if (any(is.na(X))) {
    X = matrix(NA, nrow=N, ncol=p)
    if (p>0) {
      for (f in 1:p) {
        X[,f] = sample(seq(-1, 1, length.out=N), N, replace=F)
      }
    }
  }
  
  if (p>0)
    X = scale(X)
  
  # U: image-level precision covariates
  if (any(is.na(U))) {
    U = matrix(NA, nrow=N, ncol=q)
    if (q>0) {
      for (f in 1:q) {
        U[,f] = sample(seq(-1, 1, length.out=N), N, replace=F)
      }
    }
  }
  
  if (q>0)
    U = scale(U)
  
  
  ## true underlying process ##
  
  # theta: cutoffs for Y
  theta = theta_fn(theta_tilde)
  
  # distribution of each Y
  log_pY = sapply(1:L, 
                  function(l)
                    prior('Y',
                          list(N=N, X=X), 
                          list(Y=rep(l, N), beta0=beta0, beta=beta, theta_tilde=theta_tilde)))
  pY = exp(log_pY)
  
  # sample Y: true ordinal value
  Y = apply(pY, 1, rcat, n=1)
  
  
  ## manual annotation process ##
  
  # phi: cutoffs for Z
  phi = phi_fn(phi_tilde)
  
  # randomly assign annotators to each image
  a = rep(NA, R)
  a[Z_labelled] = rcat(sum(Z_labelled), rep(1/L, A))
  
  # annotation means (nu_tilde adjusted for bias)
  nu = sapply(1:A, \(.) rnorm(L, mean=nu_tilde, sd=sqrt(nu_var))) |> t()
  
  # likelihood for Z
  log_fZ = sapply(1:L, 
                  function(l)
                    likel(list(Z=rep(l, R),  a=a, r=r, R=R, 
                               Z_labelled=Z_labelled, seq_ID=seq_ID), 
                          list(Y=Y, nu=nu, phi_tilde=phi_tilde),
                          eval_C=F, compilation='image')$log_f)
  fZ = exp(log_fZ)
  
  # sample Z: manually annotated ordinal response
  Z = rep(NA, R)
  Z[Z_labelled] = apply(fZ[Z_labelled,], 1, rcat, n=1)
  
  
  ## compositional response ##
  
  # concentrations for C
  conc_prec = concentration_fn(list(N=N, U=U, seq_ID=seq_ID),
                               list(Y=Y, alpha=alpha, omega0=omega0, omega=omega))
  concentrations_tilde = fix_zeros(conc_prec$concentration, zeta)  # fix zeros in concentration 
  
  # sample C: compositional AI predictions
  C = rdirichlet(R, concentrations_tilde)
  C_tilde = fix_zeros(C, zeta=zeta)  # fix zeros in C
  
  # whether each individual was rounded
  unrounded_indiv = rep(NA, N)
  if (any(C_labelled)) {
    unrounded_img = apply(C == C_tilde | is.na(C), 1, any)
    unrounded_indiv = sapply(1:N, function(i) all(unrounded_img[seq_ID==i]))
  }
    
  
  # observable data #
  
  N_is = N - N_oos
  
  data = list(Z=Z[1:sum(r[1:N_is])], C=C[1:sum(r[1:N_is]),], C_tilde=C_tilde[1:sum(r[1:N_is]),], 
              X=X[1:N_is,], U=U[1:sum(r[1:N_is]),], a=a[1:sum(r[1:N_is])], 
              L=L, N=N_is, r=r[1:N_is], seq_ID=seq_ID[1:sum(r[1:N_is])], R=sum(r[1:N_is]), p=p, q=q, A=A,
              Z_labelled=Z_labelled[1:sum(r[1:N_is])], C_labelled=C_labelled[1:sum(r[1:N_is])])
    
  # unobservable useful information #
  
  extra = list(Y=Y[1:N_is], 
               beta0=beta0, beta=beta, theta_tilde=theta_tilde, theta=theta, 
               nu=nu, nu_tilde=nu_tilde, phi_tilde=phi_tilde, phi=phi,
               alpha=alpha, omega0=omega0, omega=omega, 
               seed=seed)
  
  # out-of-sample #
  
  data_oos = NA
  extra_oos = NA
  if (N_oos>0) {
    data_oos = list(Z=Z[(sum(r[1:N_is])+1):R], C=C[(sum(r[1:N_is])+1):R,], C_tilde=C_tilde[(sum(r[1:N_is])+1):R,], 
                    X=X[(N_is+1):N,], U=U[(sum(r[1:N_is])+1):R,], a=a[(sum(r[1:N_is])+1):R], 
                    L=L, N=N_oos, r=r[(N_is+1):N], seq_ID=seq_ID[(sum(r[1:N_is])+1):R], R=sum(r[(N_is+1):N]), p=p, q=q, A=A,
                    Z_labelled=Z_labelled[(sum(r[1:N_is])+1):R], C_labelled=C_labelled[(sum(r[1:N_is])+1):R])
    
    extra_oos = list(Y=Y[(N_is+1):N], 
                     beta0=beta0, beta=beta, theta_tilde=theta_tilde, theta=theta,
                     nu=nu, nu_tilde=nu_tilde, phi_tilde=phi_tilde, phi=phi,
                     alpha=alpha, omega0=omega0, omega=omega, 
                  seed=seed)  
  }
  
  
  ## return ##
  
  return(list(data=data, extra=extra, data_oos=data_oos, extra_oos=extra_oos))
}


# convert data to remove ordinal/compositional responses or convert to max value/linear format
convert_data = function(data, extra=NA,  # data and extra to convert
                        new_Z_labelled=data$Z_labelled,  # new ordinal
                        new_C_labelled=data$C_labelled,  # new compositional
                        to_max=F,  # whether to convert to max value
                        to_linear=F,  # whether to convert to linear value
                        C_threshold=0.9 # threshold for which max values are acceptable
) {
  
  # fix new ordinal and compositional if necessary #
  
  if (length(new_Z_labelled)==1)
    new_Z_labelled = rep(new_Z_labelled, data$R)
  if (length(new_C_labelled)==1)
    new_C_labelled = rep(new_C_labelled, data$R)
  
  if (any(new_Z_labelled & is.na(data$Z) | 
          new_C_labelled & apply(is.na(data$C), 1, any))) {
    stop('Requested data is does not exist in the given dataset!')
  }
  
  # convert to max value/linear format #
  
  if (to_max | to_linear) {
    
    data$C[!new_C_labelled,] = NA  # remove hidden compositional images
    
    # which value is max in each row of C
    max_k = rep(NA, data$R)  
    max_k[new_C_labelled] = apply(data$C[new_C_labelled,], 1, which.max)
    
    # whether the max value meets the threshold
    max_valid = rep(FALSE, data$R)  
    max_valid[new_C_labelled] = apply(data$C[new_C_labelled,], 1, max) > C_threshold
    
    M = rep(NA, data$N)  # new "Z" value (depending on format)
    
    data$img_used = rep(F, data$R)
    
    # for each sequence...
    for (i in 1:data$N) {
      M_i = NA  # current M value
      seq_ID_i = data$seq_ID == i  # images for current sequence
      
      # convert to max value
      if (to_max) {
        max_i = max_k[seq_ID_i & max_valid]  # valid max values
        if (length(max_i) > 0) {  # if more than one, get median 
          M_i = median(max_i)
          
          # randomly round if necessary
          if ((M_i %% 1) != 0)
            M_i = ifelse(runif(1)<0.5, floor(M_i), ceiling(M_i))
          
          data$img_used[seq_ID_i & max_valid] = T
        }
      }
      
      # convert to linear value
      if (to_linear) {
        # use mean of ordinal responses if available
        Z_labelled_i = data$Z_labelled[seq_ID_i] & new_Z_labelled[seq_ID_i]  
        if (any(Z_labelled_i)) {
          M_i = mean(data$Z[seq_ID_i][Z_labelled_i])
          
          data$img_used[seq_ID_i][Z_labelled_i] = T
        } else {
          # otherwise use mean of compositional responses
          C_labelled_i = data$C_labelled[seq_ID_i] & new_C_labelled[seq_ID_i] & max_valid[seq_ID_i]
          if (any(C_labelled_i)) {
            M_i = sum((1:data$L) * colMeans(matrix(matrix(data$C[seq_ID_i,], ncol=data$L)[C_labelled_i,], ncol=data$L)))
          
            data$img_used[seq_ID_i][C_labelled_i] = T
          }
        }
        
      }
      
      M[i] = M_i  # update new "Z" value
    }
    
    # update data #
    
    data$Z_labelled = !is.na(M)
    data$C_labelled = rep(F, data$N)
    data$seq_ID = 1:data$N
    
    data$Z = M
    data$C = matrix(NA, nrow=data$N, ncol=data$L)
    data$C_tilde = data$C
    data$unrounded_indiv = rep(NA, data$N)
    
    data$r = rep(1, data$N)
    data$r[!data$Z_labelled] = 0
    data$R = sum(data$r)
    
    data$a = rep(1, data$N)
    data$a[!data$Z_labelled] = NA
    data$A = 1
    
    data$U = matrix(NA, nrow=data$N, ncol=data$s)
    
    
  } else {
  
    data$Z[!new_Z_labelled] = NA
    data$C[!new_C_labelled,] = NA
    data$C_tilde[!new_C_labelled,] = NA
    
    data$a[!new_Z_labelled] = NA
    data$Z_labelled = new_Z_labelled
    data$C_labelled = new_C_labelled
    
    data$img_used = data$Z_labelled | data$C_labelled
  }
  
  
  # remove missing data #
  
  valid_img = data$Z_labelled | data$C_labelled
  valid_indiv = rep(T, data$N)
  for (i in 1:data$N) {
    valid_img_i = valid_img[data$seq_ID==i]
    valid_indiv[i] = any(valid_img_i)
    data$r[i] = sum(valid_img_i)
  }
  data$seq_ID[!valid_img] = NA
  for (i in (1:data$N)[valid_indiv]) {
    seq_ID_i = data$seq_ID==i
    prev_i = data$seq_ID[data$seq_ID %in% 0:(i-1)]
    new_i = 1
    if (length(prev_i)>0)
      new_i = max(prev_i, na.rm=T) + 1
    data$seq_ID[seq_ID_i] = new_i
  }
  
  data$Z_labelled = data$Z_labelled[valid_img]
  data$C_labelled = data$C_labelled[valid_img]
  
  data$N = sum(valid_indiv)
  data$r = data$r[valid_indiv]
  data$R = sum(data$r)
  data$seq_ID = data$seq_ID[valid_img]
  
  data$Z = data$Z[valid_img]
  data$C = data$C[valid_img,]
  data$C_tilde = data$C_tilde[valid_img,]
  
  data$a = data$a[valid_img]
  
  data$X = as.matrix(as.matrix(data$X)[valid_indiv,], ncol=data$p)
  data$U = as.matrix(as.matrix(data$U)[valid_img,], ncol=data$q)
  
  
  unrounded_indiv = rep(NA, data$N)
  if (any(data$C_labelled)) {
    unrounded_img = apply(data$C == data$C_tilde | is.na(data$C), 1, any)
    unrounded_indiv = sapply(1:data$N, function(i) all(unrounded_img[data$seq_ID==i]))
  }
  data$unrounded_indiv = unrounded_indiv
  
  if (class(extra)=='list') {
    extra$Y = extra$Y[valid_indiv]
    extra$pY = extra$pY[valid_indiv,]
    extra$precision = extra$precision[valid_img]
    extra$concentration = extra$concentration[valid_img,]
  }
  
  # return #
  
  return(list(data=data, extra=extra))
}


## MCMC functions ##

# LIKELIHOOD #

# logged
likel = function(data, params,  # inputs
                 compilation='sum',  # how to combine likels (sum/cateogry/individual)
                 eval_Z=TRUE,  # whether to compute f_Z
                 eval_C=TRUE,  # whether to compute f_C
                 log_fZ=rep(NA, data$R), log_fC=rep(NA, data$R),  # saved previous values
                 remove_rounded=FALSE,  # whether to remove rounded compositional values
                 fallible=TRUE  # whether the annotators are fallible
) {
  
  # log_fZ=rep(NA, data$R); log_fC=rep(NA, data$R)
  
  # data
  a = data$a
  seq_ID = data$seq_ID  
  
  # params
  Y = params$Y
  Y_img = Y[data$seq_ID]
  
  ## ordinal ##
  
  if (eval_Z) {

    Z = data$Z
    A = data$A
    
    valid_Z = is.na(log_fZ) & data$Z_labelled
    
    if (fallible) {  # (default) more flexible annotators
      
      nu = params$nu  
      phi_tilde = params$phi_tilde  
      
      if (sum(valid_Z)>0) {
        
        # phi: cutoffs for ordinal responses
        phi = phi_fn(phi_tilde)
        
        # update log likel 
        log_fZ[valid_Z] = sapply((1:data$R)[valid_Z], 
                                 function(i) {
                                   pnorm(phi[Y_img[i],Z[i]+1], mean=nu[a[i],Y_img[i]]) - 
                                     pnorm(phi[Y_img[i],Z[i]], mean=nu[a[i],Y_img[i]])
                                 }) |> log()
        
      }
      
    } else {  # infallible annotators (if more than one, use proportion)
      
      for (i in (1:data$N)[valid_Z]) {
        log_fZ[i] = mean(data$Z[data$seq_ID==i] == params$Y[i]) |> 
          log()
      }

    }
  }
  
  ## compositional ##
  
  if (eval_C) {
    
    C = data$C_tilde
    U = data$U  
    
    valid_C = is.na(log_fC) & data$C_labelled 
    if (remove_rounded) {
      valid_C = valid_C & data$unrounded_indiv[data$seq_ID]
      log_fC[!data$unrounded_indiv[data$seq_ID]] = NA
    }
    
    if (sum(valid_C)>0) {
      
      alpha = params$alpha  
      omega0 = params$omega0
      omega = params$omega  
      
      # precisions and conccentrations
      conc_perp = concentration_fn(list(U=U,
                                        seq_ID=seq_ID),
                                   list(Y=Y,
                                        alpha=alpha, 
                                        omega0=omega0,
                                        omega=omega),
                                   eval=valid_C)
      
      # log likel for each image
      log_fC[valid_C] =  
        lgamma(conc_perp$precision) - 
        rowSums(lgamma(conc_perp$concentration)) + 
        rowSums((conc_perp$concentration - 1)*log(C[valid_C,]))
      
    }
  }
  
  ## joint ##
  
  # combine ordinal and compositional likelihoods
  if (eval_Z & !eval_C)
    log_f = log_fZ
  if (!eval_Z & eval_C)
    log_f = log_fC
  if (eval_Z & eval_C) {
    log_f = rep(0, data$R)
    log_f[valid_Z] = log_f[valid_Z] + log_fZ[valid_Z]
    log_f[valid_C] = log_f[valid_C] + log_fC[valid_C]
    log_f[!valid_Z & !valid_C] = NA
  }
  
  # compile into requested format
  if (compilation=='sum') {
    log_f = sum(log_f, na.rm=T)
  } else if (compilation=='category') {
    log_f = r_sum(log_f, data$r, na.rm=T)
    log_f = sapply(1:data$L, function(l) sum(log_f[Y==l], na.rm=T))
  } else if (compilation=='annotator') {
    log_f = sapply(1:data$A, function(a) sum(log_f[A==a], na.rm=T))
  } else if (compilation=='annotator category') {
    log_f2 = matrix(NA, nrow=data$A, ncol=data$L)
    for (a in 1:data$A) {
      for (l in 1:data$L) {
        log_f2[a,l] = sum(log_f[A==a & Y_img==l], na.rm=T)
      }
    }
    log_f = log_f2
  } else if (compilation=='individual') {
    log_f = r_sum(log_f, data$r, na.rm=T)
  } else if (compilation=='image') {
    
  } else {
    warning('Unknown likelihood compilation requested.')
  }
  
  
  ## return ##
  
  return(list(log_f=log_f, log_fZ=log_fZ, log_fC=log_fC))
}


# PRIOR #

# NOTE: logged; assumes cancellation when calculating acceptance probability
prior = function(.param,  # param of interest,
                 data, params,  # inputs
                 hp  # hyperparameters
) {
  
  log_p = 0  # log prior
  
  # some precalculations, if applicable
  
  if (.param=='theta_tilde') {
    default_theta = default_theta_fn(hp$prior$max_Y_prob, data$L)
    theta0 = default_theta$theta0
    theta_tilde0 = log(theta0)
  }
  
  if (.param %in% c('nu', 'nu_tilde', 'phi_tilde')) {
    default_phi = default_phi_fn(hp$prior$max_Z_prob, data$L, data$A)
    phi_tilde0 = default_phi$phi_tilde
    phi0 = default_phi$phi
  }
  
  ## Y ~ Cat(areas under curve) ##
  if (.param %in% c('Y', 'beta0', 'beta', 'theta_tilde')) {
    # Xbeta: mean for ordinal response
    Xbeta = mu_fn(data$X, params$beta0, params$beta)
    
    # theta: cutoffs for ordinal responses
    theta = theta_fn(params$theta_tilde)
    
    log_p_Z = sapply(1:data$N, 
                     function(i) {
                       pnorm(theta[params$Y[i]+1] - Xbeta[i]) - pnorm(theta[params$Y[i]] - Xbeta[i])
                     }) |> log()
    
    if (.param=='Y') {
      log_p = log_p + log_p_Z
    } else {
      log_p = log_p + sum(log_p_Z)
    }
  }
  
  ## beta0 ~ Normal(0, beta_var) ##
  if (.param=='beta0') {
    log_p = log_p + -(params$beta0)^2/(2*hp$prior$beta0_var)
  }
  
  ## beta ~ Normal(0, beta_var) ##
  if (.param=='beta') {
    log_p = log_p + sum(-(params$beta)^2/(2*hp$prior$beta_var))
  }
  
  ## theta ~ LogNormal(theta0, theta_var) ##
  if (.param=='theta_tilde') {
    log_p = log_p + sum(-(params$theta_tilde-theta_tilde0+hp$prior$theta_var/2)^2/(2*hp$prior$theta_var))
  }
  
  ## nu | nu_tilde ~ Normal(nu_tilde, nu_var) ##
  if (.param %in% c('nu', 'nu_tilde')) {
    nu_diff = apply(params$nu, 1, \(.) . - params$nu_tilde) |> t()
    log_p_nu = -(nu_diff)^2/(2*hp$prior$nu_var)
    if (.param=='nu_tilde')
      log_p_nu = colSums(log_p_nu)
    log_p = log_p + log_p_nu 
  }
  
  ## nu_tilde_l | phi_l ~ Unif(phi_{l, l-1}, phi_{l, l}) (or +/- 3 for edges) ##
  if (.param %in% c('nu_tilde', 'phi_tilde')) {
    phi = phi_fn(params$phi_tilde)
    nu_tilde_bds = nu_tilde_bds_fn(phi)
    
    log_p_nu_tilde = rep(NA, data$L)
    for (l in 1:data$L) {
      log_p_nu_tilde[l] = log(nu_tilde_bds[l,1] <= params$nu_tilde[l] & params$nu_tilde[l] <= nu_tilde_bds[l,2])
    }

    if (.param=='phi_tilde')
      log_p_nu_tilde = sum(log_p_nu_tilde)

    log_p = log_p + log_p_nu_tilde
  }
  
  ## phi ~ LogNormal(phi0, phi_var) ##
  if (.param=='phi_tilde') {
    log_p = log_p + rowSums(-(params$phi_tilde-phi_tilde0+hp$prior$phi_var/2)^2/(2*hp$prior$phi_var))
  }
  
  ## alpha ~ Dirichlet(alpha_prec*alpha0) ##
  if (.param=='alpha') {
    alpha0 = hp$prior$mu_alpha*diag(data$L)+(1-hp$prior$mu_alpha)/(data$L-1)*(1-diag(data$L))
    
    log_p = log_p + rowSums((hp$prior$alpha_precision*alpha0-1)*(log(params$alpha)))
  }
  
  ## omega0 ~ Normal(omega0_mn, omega0_var) ##
  if (.param=='omega0') {
    log_p = log_p + -(params$omega0 - hp$prior$omega0_mn)^2/(2*hp$prior$omega_var)
  }
  
  ## omega ~ Normal(0, omega_var) ##
  if (.param=='omega') {
    log_p = log_p + sum(-(params$omega)^2/(2*hp$prior$omega_var))
  }
  
  ## return ##
  
  return(log_p)
  
}


# PROPOSALS #

Y_Gibbs_fn = function(data, params, 
                           log_fZ=rep(NA, data$R), log_fC=rep(NA, data$R),  # saved previous values
                           fallible=T,  # whether the annotators are fallible
                           max_value=700  # threshold for logged numerator values
                           ) {
  
  # likelihoods for Z, C, Y for each sequence and possible category
  Z_likel = matrix(NA, nrow=data$N, ncol=data$L)
  C_likel = matrix(NA, nrow=data$N, ncol=data$L)
  Y_likel = matrix(NA, nrow=data$N, ncol=data$L)
  
  # likelihoods for Z, C for each image and possible category (for future computations)
  Z_likel_img = matrix(NA, nrow=data$R, ncol=data$L)
  C_likel_img = matrix(NA, nrow=data$R, ncol=data$L)
  
  # helper theoretical parameters
  params2 = params
  
  # get likelihoods for each category
  for (l in 1:data$L) {
    
    # update params
    params2$Y = rep(l, data$N)
    
    # likelihoods of Z and C #
    
    # use saved likelihoods for Y that match current params
    Y_match = params$Y == l
    Y_img_match = Y_match[data$seq_ID]
    Y_img_match[is.na(Y_img_match)] = F
    
    log_fZ_l = rep(NA, data$R)
    log_fZ_l[Y_img_match] = log_fZ[Y_img_match]
    
    log_fC_l = rep(NA, data$R)
    log_fC_l[Y_img_match] = log_fC[Y_img_match]
    
    # compute likelihoods
    Z_likel_l = likel(data, params2, eval_C=F, compilation='individual',
                      log_fZ=log_fZ_l, log_fC=log_fC_l, fallible=fallible)
    C_likel_l = likel(data, params2, eval_Z=F, compilation='individual',
                      log_fZ=log_fZ_l, log_fC=log_fC_l, fallible=fallible)
    
    # update matrices
    Z_likel[,l] = Z_likel_l$log_f
    C_likel[,l] = C_likel_l$log_f
    
    Z_likel_img[,l] = Z_likel_l$log_fZ
    C_likel_img[,l] = C_likel_l$log_fC
    
    # likelihood of Y #
    
    Y_likel[,l] = prior('Y', data, params2)
    
  }
  
  # compute distribution of Y | Z, C, beta, ...
  num = C_likel + Z_likel + Y_likel  # logged numerator
  max_value = 700  # computational fix
  num_bad_vals = rowSums(num>max_value)
  if (max(num_bad_vals)>0) {
    warning_message = paste(sum(num_bad_vals==1), ' rows with single value and ',
                            sum(num_bad_vals>1), ' rows with more than one value above ',
                            max_value, ' thresholded.', 
                            sep='')
    warning(warning_message)
  }
  num[num>max_value] = max_value
  exp_num = exp(num)
  denom = log(rowSums(exp_num))  # logged denominator to normalize
  Y_prob = exp(num - denom)  # distribution
  
  return(list(Y_prob=Y_prob,
              Z_likel_img=Z_likel_img,
              C_likel_img=C_likel_img))
}

# NOTE: logged; assumes cancellation when calculating acceptance probability
proposals = function(.param,  # param of interest
                     data, params,  # inputs
                     hp,  # hyperparameters
                     log_fZ=rep(NA, data$R), log_fC=rep(NA, data$R),  # saved previous values
                     fallible=T,  # whether the annotators are fallible
                     seed=round(runif(1, min=1, max=1e9))  # random seed 
) {
  
  # set seed for reproducibility
  set.seed(seed)
  
  params_new = params  # proposed parameters
  log_prop_diff = NA  # log difference in proposal probabilities (old to new - new to old)
  
  ## new Z using Gibbs ##
  if (.param=='Y') {
    
    Y_Gibbs_results = Y_Gibbs_fn(data, params, 
                                           log_fZ=log_fZ, log_fC=log_fC, fallible=fallible)
    
    # sample new Y
    params_new$Y = apply(Y_Gibbs_results$Y_prob, 1, rcat, n=1)
    
    # update saved likelihoods
    for (i in 1:data$R) {
      Y_i = params_new$Y[data$seq_ID][i]
      log_fZ[i] = Y_Gibbs_results$Z_likel_img[i,Y_i]
      log_fC[i] = Y_Gibbs_results$C_likel_img[i,Y_i]
    }
  
  }
  
  ## new beta0 ~ Normal(prev beta0, hp$prop$beta0) ##
  if (.param=='beta0') {
    params_new$beta0 = rnorm(1, mean=params$beta0, sd=hp$prop$beta0)
    log_prop_diff = 0  # symmetric
  }
  
  ## new beta ~ Normal(prev beta, hp$prop$beta) ##
  if (.param=='beta') {
    params_new$beta = rnorm(data$p, mean=params$beta, sd=hp$prop$beta)
    log_prop_diff = 0  # symmetric
  }
  
  ## new theta_tilde ~ Normal(prev theta_tilde, hp$prop$theta_tilde) ##
  if (.param=='theta_tilde') {
    params_new$theta_tilde = rnorm(data$L-2, mean=params$theta_tilde-hp$prop$theta_tilde^2/2, sd=hp$prop$theta_tilde)
    log_prop_diff = 0  # symmetric
  }
  
  ## new nu ~ Normal(prev nu, hp$prop$nu) ##
  if (.param=='nu') {
    params_new$nu = rnorm(data$A*data$L, mean=c(params$nu), sd=hp$prop$nu) |> 
      matrix(ncol=data$L)
    log_prop_diff = 0  # symmetric
  }
  
  ## new nu_tilde_l ~ TruncNormal(prev nu_tilde, hp$prop$nu_tilde, within bounds of phi) ##
  if (.param=='nu_tilde') {
    phi = phi_fn(params$phi_tilde)
    nu_tilde_bds = nu_tilde_bds_fn(phi)
    
    log_prop_diff = rep(NA, data$L)
    
    for (l in 1:data$L) {
      params_new$nu_tilde[l] = rtrunc(1, 'norm', 
                                   a=nu_tilde_bds[l,1], b=nu_tilde_bds[l,2],
                                   mean=params$nu_tilde[l], sd=hp$prop$nu_tilde[l])
      
      log_prop_diff[l] = 
        dtrunc(params_new$nu_tilde[l], 'norm', 
               a=nu_tilde_bds[l,1], b=nu_tilde_bds[l,2],
               mean=params$nu_tilde[l], sd=hp$prop$nu_tilde[l]) -
        dtrunc(params$nu_tilde[l], 'norm', 
               a=nu_tilde_bds[l,1], b=nu_tilde_bds[l,2],
               mean=params_new$nu_tilde[l], sd=hp$prop$nu_tilde[l])
    }
  }
  
  ## new phi_tilde ~ Normal(prev phi_tilde, hp$prop$phi_tilde) ##
  if (.param=='phi_tilde') {
    phi_tilde_mn = apply(params$phi_tilde, 2, \(.) .-hp$prop$phi_tilde^2/2) 
    params_new$phi_tilde = rnorm(data$L*(data$L-2), mean=c(phi_tilde_mn), sd=rep(hp$prop$phi_tilde, data$L-2)) |> 
      matrix(nrow=data$L)
    log_prop_diff = 0  # symmetric
  }
  
  ## new alpha ~ Dirichlet(hp$prop$alpha * prev alpha) ##
  if (.param=='alpha') {
    params_new$alpha = rdirichlet(data$L, hp$prop$alpha * params$alpha)
    
    log_prop_diff = 
      sapply(1:data$L, 
             function(l)
               (-sum(lgamma(hp$prop$alpha[l]*params$alpha[l,])) + sum((hp$prop$alpha[l]*params$alpha[l,]-1)*log(params_new$alpha[l,]))) -
               (-sum(lgamma(hp$prop$alpha[l]*params_new$alpha[l,])) + sum((hp$prop$alpha[l]*params_new$alpha[l,]-1)*log(params$alpha[l,])))
      )
  }
  
  ## new omega0 ~ Normal(prev omega0, hp$prop$omega0) ##
  if (.param=='omega0') {
    params_new$omega0 = rnorm(1, mean=params$omega0, sd=hp$prop$omega0)
    log_prop_diff = 0  # symmetric
  }
  
  ## new omega ~ Normal(prev omega, hp$prop$omega) ##
  if (.param=='omega') {
    params_new$omega = rnorm(data$q, mean=params$omega, sd=hp$prop$omega)
    log_prop_diff = 0  # symmetric
  }
  
  ## return ##
  
  return(list(params_new=params_new, log_prop_diff=log_prop_diff,
              log_fZ=log_fZ, log_fC=log_fC))
  
}


# UPDATE PARAMETERS #

# make proposals, update parameters and acceptances
update_params = function(.param,  # parameter of interest
                         data, params,  # inputs
                         hp,  # hyperparameters
                         log_fZ=rep(NA, data$R), log_fC=rep(NA, data$R),  # saved previous values
                         fallible=T,  # whether the annotators are fallible
                         label='',  # label for error message
                         seed=round(runif(1, min=1, max=1e9))  # random seed 
) {
  
  ## proposal ##
  
  prop = proposals(.param, data, params, hp, 
                   log_fZ=log_fZ, log_fC=log_fC, fallible=fallible,
                   seed=seed)
  # proposed parameters
  params_new = prop$params_new  
  # log difference in proposal probabilities (old to new - new to old)
  log_prop_diff = prop$log_prop_diff  
  
  if (.param=='Y') {
    
    # update likelihood
    log_fZ = prop$log_fZ
    log_fC = prop$log_fC
    
    
    ## return ##
    
    return(list(params=params_new,
                log_fZ=log_fZ, log_fC=log_fC))
    
  } else {
    
    ## prior ##
    
    # log difference in prior probabilities (new - old)
    log_prior_diff = 
      prior(.param, data, params_new, hp) -
      prior(.param, data, params, hp)
    
    ## likelihood ##
    
    log_likel_diff = 0
    eval_Z = F; eval_C = F
    if (.param %in% c('nu', 'phi_tilde', 'alpha', 'omega0', 'omega')) {
      
      # whether to evaluate ordinal, compositional, or both likelihoods
      eval_Z = .param %in% c('nu', 'phi_tilde')
      eval_C = .param %in% c('alpha', 'omega0', 'omega')
      
      # whether to sum likelihood or report for each individual
      compilation = 'sum'
      if (.param %in% c('phi_tilde', 'alpha'))
        compilation = 'category'
      if (.param %in% c())
        compilation = 'annotator'
      if (.param %in% c('nu'))
        compilation = 'annotator category'
      if (.param %in% c())
        compilation = 'individual'
      
      # whether to exclude rounded values of C
      remove_rounded = F #.param %in% c('omega0', 'omega')
      
      # log difference in likelihoods (new - old)
      log_likel_new = likel(data, params_new, compilation=compilation, 
                            eval_Z=eval_Z, eval_C=eval_C,
                            remove_rounded=remove_rounded,
                            fallible=fallible)
      log_likel_old = likel(data, params, compilation=compilation, 
                            eval_Z=eval_Z, eval_C=eval_C,
                            log_fZ=log_fZ, log_fC=log_fC,
                            remove_rounded=remove_rounded,
                            fallible=fallible)
      
      log_likel_diff = 
        log_likel_new$log_f -
        log_likel_old$log_f
  
    }
    
    ## acceptance ##
    
    # calculate ratio used to determine acceptance
    log_R = log_likel_diff + log_prior_diff - log_prop_diff
    R = as.vector(exp(log_R))  
    R = sapply(R, function(x) min(x, 1))
    n_to_update = length(R)  # number of acceptance checks
    
    # random uniform sample(s)
    set.seed(seed*1e2)
    S = runif(n_to_update)  
    
    # check whether the proposal(s) are accepted
    accepted = S < R  
    
    # update if accepted 
    if (any(is.na(accepted)) | any(is.null(accepted))) {
      error_log = list(.param=.param,
                       params=params,
                       params_new=params_new,
                       data=data,
                       hp=hp,
                       seed=seed)
      filename = paste('error_log_', label, '.RData', sep='')
      save(error_log, file=filename)
      stop(paste('Something went wrong in checking acceptance. Current settings saved in ',
                 paste(getwd(), filename, sep='/'), sep=''))
    }
    
    params[[.param]][accepted] = params_new[[.param]][accepted]
    
    
    # update likelihoods
    if (.param %in% c('omega0', 'omega')) {
      if (accepted) {
        if (eval_Z)
          log_fZ = log_likel_new$log_fZ
        if (eval_C)
          log_fC = log_likel_new$log_fC
      }
    } else if (.param == 'nu') {
      accepted = matrix(accepted, nrow=data$A, ncol=data$L)
      Y_img = params$Y[data$seq_ID]
      for (i in (1:data$R)[data$Z_labelled]) {
        if (accepted[data$a[i], Y_img[i]])
          log_fZ[i] = log_likel_new$log_fZ[i]
      }
    } else if (.param %in% c('phi_tilde', 'alpha')) {
      accepted_l = (1:data$L)[accepted]
      Y_l = params_new$Y %in% accepted_l
      if (eval_Z)
        log_fZ[Y_l[data$seq_ID]] = log_likel_new$log_fZ[Y_l[data$seq_ID]]
      if (eval_C)
        log_fC[Y_l[data$seq_ID]] = log_likel_new$log_fC[Y_l[data$seq_ID]]
    }
    
    ## return ##
    
    return(list(params=params, accepted=accepted,
                log_fZ=log_fZ, log_fC=log_fC))
    
  } 
  
}


# MCMC (METROPOLIS-HASTINGS) #

MH_alg = function(n_iter, # number of iterations
                  data,  # named list of observed data
                  params0=NA,  # named list of initial parameters 
                  sampling_params=NA,  # parameters to be sampled
                  fixed_params=NA,  # parameters to be fixed
                  thin=1,  # thinning
                  prior0=NA,  # names list of hyperparameters to priors
                  prop0=NA,  # named list of initial hyperparameters for proposals
                  fallible=T,  # whether the annotators are fallible
                  adaptive_sampling=FALSE,  # adaptive MCMC
                  cluster_size=100,  # size of cluster for adaptive sampling
                  verbose=TRUE,  # whether to print notifications
                  progressbar=TRUE,  # whether to include a progressbar
                  label=round(proc.time()[3]),  # label for error message
                  seed=round(runif(1, min=1, max=1e9))  # random seed
) {
  
  ## initial parameters ##
  
  params0 = default_params0_fn(data, params0)
  
  ## sampling params ##
  
  sampling_params = default_sampling_params_fn(data, 
                                               sampling_params=sampling_params, 
                                               fixed_params=fixed_params,
                                               fallible=fallible)
  
  ## structures to hold results ##
  
  # generate names for parameters to be sampled
  param_names = default_param_names_fn(data)
  
  # generate empty structures for samples and acceptance rates
  samples = list()
  acceptance_rates = list()
  for (.param in sampling_params) {
    sampling_pn = param_names$sampling_param_names[[.param]]
    samples[[.param]] = matrix(NA, nrow=n_iter, ncol=length(sampling_pn))
    colnames(samples[[.param]]) = sampling_pn
    if (.param!='Y') {
      acceptance_pn = param_names$acceptance_param_names[[.param]]
      acceptance_rates[[.param]] = rep(0, length(acceptance_pn))
      names(acceptance_rates[[.param]]) = acceptance_pn
    }
  }
  
  cluster_acceptance_rates = acceptance_rates
  
  
  ## initialize ##
  
  # hyperparameters
  hp = default_hp_fn(data, prior0, prop0, sampling_params)
  
  # initial parameters
  params = params0
  
  # initial likelihoods
  log_likel = likel(data, params)
  log_fZ = log_likel$log_fZ
  log_fC = log_likel$log_fC
  
  # indicate which parameters are being sampled
  if (verbose) 
    print(paste('Sampling: ', paste(sampling_params, collapse=', '), sep=''))
  
  
  if (n_iter > 0) { 
    # setup progressbar
    pb = NA
    if (progressbar)
      pb = txtProgressBar(max=thin*n_iter, style=3)
    
    ## generate samples ##
    
    for (s in 1:(thin*n_iter)) {
      
      seed_s = seed + 1e3*(s-1)
      
      # iterate through parameters
      for (.param in sampling_params) {
        # update parameters 
        update = update_params(.param, data, params, hp, 
                               log_fZ=log_fZ, log_fC=log_fC, fallible=fallible,
                               label=label, seed=seed_s)
        params = update$params
        log_fZ = update$log_fZ
        log_fC = update$log_fC
        
        # update samples and acceptance rates
        accepted = update$accepted
        
        if (s%%thin==0) {
          if (.param!='Y') {
            ar = accepted / n_iter  # acceptance rate for current param
            acceptance_rates[[.param]] = acceptance_rates[[.param]] + ar
          } 
          
          s_t = s/thin
          samples[[.param]][s_t,] = as.vector(params[[.param]])
        }
        
        # adaptive sampling
        if (adaptive_sampling & .param!='Y') {
          cluster_ar = accepted / cluster_size
          cluster_acceptance_rates[[.param]] = cluster_acceptance_rates[[.param]] + cluster_ar
          
          if (s%%cluster_size==0) {
            low_ar = cluster_acceptance_rates[[.param]] < 0.2
            high_ar = cluster_acceptance_rates[[.param]] > 0.4
            
            sd_adjust = s^-0.3
            if (.param=='alpha')
              sd_adjust = -sd_adjust
            
            hp$prop[[.param]][low_ar] = hp$prop[[.param]][low_ar] * (1-sd_adjust)
            hp$prop[[.param]][high_ar] = hp$prop[[.param]][high_ar] * (1+sd_adjust)
            
            if (s!=thin*n_iter) {
              cluster_acceptance_rates[[.param]] = rep(0, length(cluster_acceptance_rates[[.param]]))
              cluster_acceptance_rates[[.param]][is.na(acceptance_rates[[.param]])] = NA 
            }
          }
        }
        
        # update seed
        seed_s = seed_s + 1e2
      }
      
      # update progressbar
      if (progressbar)
        setTxtProgressBar(pb, s)
    }
    
    
    # close progressbar
    if (progressbar)
      close(pb)
    
  } else {
    if (verbose)
      print('No iterations requested -- returning without sampling.')
  }
  
  ## return ##
  
  return(list(samples=samples,
              params=params,
              acceptance_rates=acceptance_rates,
              cluster_acceptance_rates=cluster_acceptance_rates,
              sampling_params=sampling_params,
              n_iter=n_iter,
              thin=thin,
              params0=params0,
              hp=hp,
              adaptive_sampling=adaptive_sampling,
              cluster_size=cluster_size,
              fallible=fallible,
              log_fZ=log_fZ, log_fC=log_fC,
              seed=seed))
  
}



## LINEAR MODEL ##

linear_model = function(data, n_iter, n_burn=0, seed=round(runif(1, min=1, max=1e9))) {
  
  # set seed for reproducibility
  set.seed(seed)
  
  model_str = 'model {
        # likelihood
        for(i in 1:N){
          Z[i] ~ dnorm(beta0 + inprod(X[i,],beta[]), taue)
        }
        # priors
        beta0 ~ dnorm(0, 1/10)
        for(f in 1:p){
          beta[f] ~ dnorm(0, taub*taue)
        }
        taue ~ dgamma(0.1, 0.1)
        taub ~ dgamma(0.1, 0.1)
      }'
  
  linfit = jags.model(textConnection(model_str),
                      data = list(Z=data$Z, X=data$X, N=data$N, p=data$p),
                      n.chains = 1)
  update(linfit, n_burn)
  coda_results = coda.samples(linfit, 
                              variable.names=c('beta0', 'beta'),
                              n.iter = n_iter)
  
  beta0_results = as.data.frame(coda_results[[1]])[,(data$p+1)] |> matrix(ncol=1)
  beta_results = as.data.frame(coda_results[[1]])[,-(data$p+1)] |> as.matrix()
  colnames(beta0_results) = default_param_names_fn(data)$sampling_param_names$beta0
  colnames(beta_results) = default_param_names_fn(data)$sampling_param_names$beta
  
  results = list(samples=list(beta0=beta0_results,
                              beta=beta_results),
                 n_iter=n_iter,
                 seed=seed)
  
  return(results)
}


#######################
# PERFORMANCE METRICS #
#######################

## calculate RPS given the true and predicted CDF's ##
rps_fn = function(cdf_true, # true CDF; n x L array
                  cdf_pred,  # predicted CDF; n x L x n_iter
                  to_mean=T  # whether to report mean (vs for each individual)
                  ) {
  
  L = ncol(cdf_true)
  n_iter = dim(cdf_pred)[3]
  
  cdf_pred_mn = apply(cdf_pred, 1:2, mean)  # mean over iterations (n x L)
  
  rps = rowSums((cdf_pred_mn[,-L] - cdf_true[,-L])^2) / (L-1)  # RPS for each individual
  
  if (to_mean)
    rps = mean(rps)
  
  return(rps)
}

## convert PMF to CDF ##
pmf_to_cdf = function(pmf) {
  cdf = pmf
  for (col in 2:ncol(pmf)) {
    cdf[,col] = cdf[,col-1] + cdf[,col]
  }
  return(cdf)
}

## convert categorical to CDF ##
categorical_to_cdf = function(cats, L) {
  cdf = matrix(0, nrow=length(cats), ncol=L)
  for (l in 1:L) {
    is_l = cats==l
    cdf[is_l,l:L] = 1
  }
  return(cdf)
}

## get parameters from iteration s ##
slice_sample = function(samples, s) {
  params = sapply(samples, \(param) param[s,])
  L = length(params$theta_tilde) + 2
  if (!is.null(params$nu))
    params$nu = matrix(params$nu, ncol=L)
  if (!is.null(params$phi_tilde))
    params$phi_tilde = matrix(params$phi_tilde, nrow=L)
  if (!is.null(params$alpha))
    params$alpha = matrix(params$alpha, nrow=L)
  return(params)
}

## get RPS for P(Y | Z, C) ##
Y_rps_is_fn = function(data, extra, results) {
  
  cdf_true = categorical_to_cdf(extra$Y, data$L)
  
  cdf_pred = array(NA, dim=c(data$N, data$L, results$n_iter))
  
  pb = txtProgressBar(max=results$n_iter, style=3)
  for (s in 1:results$n_iter) {
    cdf_pred[,,s] = categorical_to_cdf(results$samples$Y[s,], data$L)
    
    setTxtProgressBar(pb, s)
  }
  close(pb)
  
  rps = rps_fn(cdf_true, cdf_pred)
  
  return(rps)
}

## get RPS for P(Y_0 | beta0, beta, theta) ##
Y_rps_oos_fn = function(data, extra, results) {
  
  cdf_true = categorical_to_cdf(extra$Y, data$L)
  
  cdf_pred = array(NA, dim=c(data$N, data$L, results$n_iter))
  
  pb = txtProgressBar(max=results$n_iter, style=3)
  for (s in 1:results$n_iter) {
    params = slice_sample(results$samples, s)
    
    # distribution of each Y
    log_pY = sapply(1:data$L, 
                         function(l)
                           prior('Y',
                                 data, 
                                 list(Y=rep(l, data$N), beta0=params$beta0, beta=params$beta, theta_tilde=params$theta_tilde)))
    pY = exp(log_pY)
    
    
    # sample Y: true ordinal value
    Y_samps = apply(pY, 1, rcat, n=1)
    
    cdf_pred[,,s] = categorical_to_cdf(Y_samps, data$L)
    
    setTxtProgressBar(pb, s)
  }
  close(pb)

  rps = rps_fn(cdf_true, cdf_pred)
  
  return(rps)
}




## plotting ##

# SAMPLING DIAGNOSTICS #

# pretty print acceptance rates
print_acceptance_rates = function(acceptance_rates) {
  for (.param in names(acceptance_rates)) {
    ar = round(acceptance_rates[[.param]] * 100)
    if (length(ar)==1) {
      print(paste(.param, ': ', ar, '%', sep=''))
    } else {
      ar = ar[!is.na(ar)]
      print(paste(.param, ': mean ', round(mean(ar)), 
                  '% (min ', min(ar), '%, max ', max(ar), '%)', sep=''))
    }
  }
}

# calculate mean "confusion" matrices
confusion_matrix_fn = function(data, results, decimals=NA) {
  
  annotator_confusions = NA; combined_annotator_confusions = NA
  if ('nu_tilde' %in% results$sampling_params &
      'nu' %in% results$sampling_params & 
      'phi_tilde' %in% results$sampling_params) {
    annotator_confusions = array(0, 
                                 dim=c(data$L, data$L, data$A),
                                 dimnames=list(Y=1:data$L, 
                                               Z=1:data$L,
                                               annotator=1:data$A))
    for (s in 1:results$n_iter) {
      phi_s = phi_fn(matrix(results$samples$phi_tilde[s,], nrow=data$L))
      nu_s = matrix(results$samples$nu[s,], nrow=data$A) 
      
      for (Y in 1:data$L) {
        for (Z in 1:data$L) {
          annotator_confusions[Y,Z,] = annotator_confusions[Y,Z,] + 
            (pnorm(phi_s[Y,Z+1], mean=nu_s[,Y], sd=1) -
               pnorm(phi_s[Y,Z], mean=nu_s[,Y], sd=1))
        }
      }
    }
    
    annotator_confusions = annotator_confusions / results$n_iter
    
    combined_annotator_confusions = apply(annotator_confusions, 1:2, mean)
  }
  
  AI_confidence_confusions = NA
  if ('alpha' %in% results$sampling_params &
      'omega0' %in% results$sampling_params) {
    AI_confidence_confusions = results$samples$alpha |> 
      colMeans() |> 
      matrix(nrow=data$L)
  }
  
  if (!is.na(decimals)) {
    combined_annotator_confusions = round(combined_annotator_confusions, decimals)
    annotator_confusions = round(annotator_confusions, decimals)
    AI_confidence_confusions = round(AI_confidence_confusions, decimals)
  }
  
  
  return(list(combined_annotator=combined_annotator_confusions,
              annotator=annotator_confusions,
              AI_confidence=AI_confidence_confusions))
}

# helpful info for sequence: Y, MCMC prediction proportions, annotations, AI predictions
seq_info_fn = function(seq, data, extra=NA, samples=NA, decimals=NA, print=FALSE) {
  Y = NA
  if (class(extra)=='list')
    Y = extra$Y[seq]
  
  prediction = NA
  if (!any(is.na(samples))) {
    prediction = table(factor(samples$Y[,seq], levels=1:data$L))
    prediction = prediction / sum(prediction)
  }
  
  response = cbind(data$Z[data$seq_ID==seq],
                   matrix(data$C[data$seq_ID==seq,],
                          ncol=data$L))
  colnames(response) = c('annotation',
                         paste('conf', 1:data$L, sep='_'))
  
  if (!is.na(decimals)) {
    prediction = round(prediction, decimals)
    response = round(response, decimals)
  }
  
  if (print) {
    if (!is.na(Y))
      print(paste('Y: ', Y, sep=''))
    if (!any(is.na(prediction)))
      print(paste('prediction: [', paste(prediction, collapse=', '), ']'))
    print(knitr::kable(response))
  } else {
    return(list(Y=Y, 
                prediction=prediction,
                response=response))
  }
}

# wrapper function for trace and density plots
diagnostic_plot = function(plot_type, 
                           samples, data, extra=list(), 
                           params=NA, burnin=0, max_plot=3, plots_across=1,
                           main=list()) {
  
  # iterate through parameters
  plotting_params = names(samples)
  if (!all(is.na(params)))
    plotting_params = intersect(plotting_params, params)
  for (.param in plotting_params) {
    # extract samples for current parameter 
    param_samples = as.matrix(samples[[.param]])
    n_iter = nrow(param_samples)
    pn = colnames(samples[[.param]])
    
    if (burnin>=n_iter)
      stop('Burnin is greater than the number of iterations!')
    
    param_samples = as.matrix(param_samples[(burnin+1):n_iter,])
    n_iter = n_iter - burnin
    
    # extract number of parameters and which columns to plot, if applicable
    if (.param=='alpha') {
      to_plot = 1:data$L
      n_params = data$L
    } else {
      n_params = ncol(param_samples)
      
      params_valid = 1:n_params
      
      
      if(! .param %in% c('beta', 'omega'))
        n_params = min(max_plot, n_params)
      
      to_plot = sample(params_valid, n_params, replace=FALSE) |> sort()
    }
    
    # iterate through indices of parameter, if applicable
    across = 1
    down = length(to_plot)
    if (plot_type=='tracedensity') 
      across = 2
    across = across * plots_across
    down = down / plots_across
    par(mfrow=c(down, across), mar=c(5,4,4,2)+0.1, oma=c(0,0,0,0))
    for (i in to_plot) {
      # get current samples
      cur_samples = param_samples[,i]
      if (.param=='alpha')
        cur_samples = param_samples[,seq(i, data$L^2, by=data$L)]
      
      if (!is.null(main[[.param]])) {
        cur_main = main[[.param]][i]
      } else {
        cur_main = pn[i]
        if (.param %in% c('beta', 'omega')) {
          covar = NA
          if (.param=='beta' & !is.null(colnames(data$X)[i]))
            covar = colnames(data$X)[i]
          if (.param=='omega' & !is.null(colnames(data$U)[i]))
            covar = colnames(data$U)[i]
          if (!is.na(covar))
            cur_main = paste(cur_main, ' (', covar, ')', sep='')
        }
      }
      
      # traceplot
      if (plot_type=='trace' | plot_type=='tracedensity') {
        
        if (.param=='alpha') {
          
          plot(burnin:(burnin+n_iter-1), 
               xlim=c(burnin, burnin+n_iter-1),
               ylim=c(0, 1), 
               xlab='iteration', ylab='value', main=cur_main,
               type = "n")
          
          col = c('red', 'orange', 'green', 'blue', 'purple')
          for (l2 in 1:data$L) {
            lines(burnin:(burnin+n_iter-1), cur_samples[,l2], 
                  col=col[l2], lty=1)
            
            if (!is.null(extra$alpha[i,]))
              abline(h=extra$alpha[i,l2], col=col[l2], lty=2)
          }
          
        } else {
          plot(burnin:(burnin+n_iter-1), 
               xlim=c(burnin, burnin+n_iter-1),
               ylim=range(cur_samples), 
               xlab='iteration', ylab='value', main=cur_main,
               type = "n")
          
          lines(burnin:(burnin+n_iter-1), cur_samples, 
                col='red', lty=1)
          
          if (!is.null(extra[[.param]][i]))
            abline(h=extra[[.param]][i], lty=2)
        }
        
      } 
      
      # density plot
      if (plot_type=='density'| plot_type=='tracedensity') {
        if (.param=='Y') {
          
          cur_samples = cur_samples |> ordered(levels=1:data$L)
          
          col = rep('gray', data$L)
          if (!is.null(extra$Y[i])) {
            col = rep('red', data$L)
            col[extra$Y[i]] = 'blue'
          }
          
          barplot(table(cur_samples), col=col, space=0)
          
        } else if (.param=='alpha') {
          
          alpha_quants = apply(cur_samples, 2, quantile, prob=c(0.025, 0.5, 0.975))
          
          if (!any(is.null(extra$alpha))) {
            barplot(extra$alpha[i,], space=0, ylim=c(0, 1))
          } else {
            barplot(rep(0, data$L), space=0, ylim=c(0, 1))
          }
          arrows(x0=1:data$L-0.5, x1=1:data$L-0.5, y0=alpha_quants[1,], y1=alpha_quants[3,], 
                 code=3, length=0.1, angle=90)
          points(1:data$L-0.5, alpha_quants[2,], pch=16)

        } else {
          
          xlim = quantile(cur_samples, c(0.001, 0.999))
          xlim = range(xlim, extra[[.param]][i])
          
          # plot(density(cur_samples), 
          #      col='black', lty=1, lwd=1, xlim=xlim,
          #      xlab='value', ylab='', main=cur_main,
          #      yaxt='n', frame.plot=F)
          
          
          plot(density(cur_samples), 
               xlim=xlim, xlab='', ylab='', main=cur_main,
               xaxt='n', yaxt='n', frame.plot=F, type='n')
          polygon(density(cur_samples), 
               border='black', col='orchid1', lty=1, lwd=2)
          # plot(density(cur_samples), 
          #      col='black', lty=1, lwd=2, xlim=xlim,
          #      xlab='', ylab='', main=cur_main,
          #      xaxt='n', yaxt='n', frame.plot=F)
          axis(side=1, lwd=2, pos=0)
          
          if (!is.null(extra[[.param]][i]))
            abline(v=extra[[.param]][i], col='blue', lty=2)
          
        } 
      }
      
      # acf plot
      if (plot_type=='acf') {
        acf(cur_samples, main=pn[i])
      }
    }
  }
  
  par(mfrow=c(1,1))
}

# traceplots
traceplots = function(samples, data, extra=list(), params=NA, burnin=0, max_plot=3, plots_across=1, main=list()) {
  diagnostic_plot('trace', samples, data, extra, params=params, burnin=burnin, max_plot=max_plot, plots_across=plots_across, main=main)
}

# density plots 
densplots = function(samples, data, extra=list(), params=NA, burnin=0, max_plot=3, plots_across=1, main=list()) {
  diagnostic_plot('density', samples, data, extra, params=params, burnin=burnin, max_plot=max_plot, plots_across=plots_across, main=main)
}

# traceplots and density plots combined
tracedensplots = function(samples, data, extra=list(), params=NA, burnin=0, max_plot=3, plots_across=1, main=list()) {
  diagnostic_plot('tracedensity', samples, data, extra, params=params, burnin=burnin, max_plot=max_plot, plots_across=plots_across, main=main)
}

# acf plot
acfplots = function(samples, data, extra=list(), params=NA, burnin=0, max_plot=3, plots_across=1, main=list()) {
  diagnostic_plot('acf', samples, data, extra, params=params, burnin=burnin, max_plot=max_plot, plots_across=plots_across, main=main)
}



















