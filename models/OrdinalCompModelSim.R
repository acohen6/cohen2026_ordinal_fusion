
################################################################################

# Cohen et al. 2026                                                 
# Code to run simulation study

################################################################################


## imports ##

library(parallel)
library(foreach)
library(doSNOW)

## setup ##

source('models/OrdinalCompModelSetup.R')

# WARNING: make sure load variable is correct -- will rewrite files
if (!load) 
  write.table(matrix(c('setting', 'version', 'epoch', 'time'), nrow=1), here::here(save_path, 'runtime_results.csv'), row.names=F, col.names=F, sep=',')

if (!load | calc_RPS) 
  write.table(matrix(c('setting', 'version', 'epoch', 'Y_0', 'Y'), nrow=1), here::here(save_path, 'RPS_results.csv'), row.names=F, col.names=F, sep=',')


## simulations ##

n_cores = detectCores()
n_cores = min(n_epochs, n_cores - 1) 
cluster = makeCluster(n_cores)
registerDoSNOW(cluster)

pb = txtProgressBar(max=n_epochs, style=3)
progress = function(n) setTxtProgressBar(pb, n)
opts = list(progress=progress)

foreach(epoch=1:n_epochs, .options.snow=opts, .combine=rbind) %dopar% {
  
  source(sim_func_path)
  
  seed_epoch = starting_seed + epoch*1e3
  
  ## current datasets ##
  
  if (load) {
    load(here::here(save_path, sprintf('data_extra_full_%i.RData', epoch)))
    load(here::here(save_path, sprintf('max_Z_C_labelled_%i.RData', epoch)))
    
    data_is = data_extra_full$data
    extra_is = data_extra_full$extra
    
    if (!real_data) {
      data_oos = data_extra_full$data_oos
      extra_oos = data_extra_full$extra_oos
    }
  } else {
    if (!real_data) {
      # simulate data
      data_extra_full = gen_sim_data(N_is=N_is, N_oos=N_oos, covariates=covariates_all, seed=seed_epoch)
      
      data_oos = data_extra_full$data_oos
      extra_oos = data_extra_full$extra_oos
      
      # label certain parts of dataset
      
      max_Z_C_labelled_lst = list()
      # set seed for reproducibility
      set.seed(seed_epoch + 1e2)
      
      max_Z_labelled = rep(F, data_extra_full$data$R)  # current version of whether each image is annotated (if available to setting)
      for (prop_i in seq_along(prop_Z_labelled_lst)) {
        prop = prop_Z_labelled_lst[prop_i]
        
        n_to_label = ceiling(prop*data_extra_full$data$R - sum(max_Z_labelled))
        
        max_Z_labelled[sample((1:data_extra_full$data$R), n_to_label, replace=F)] = T
        
        max_C_labelled = rep(T, data_extra_full$data$R)
        
        max_Z_C_labelled_lst[[prop_i]] = list(Z=max_Z_labelled, C=max_C_labelled)
      }
      
      save(max_Z_C_labelled_lst, file=here::here(save_path, sprintf('max_Z_C_labelled_%i.RData', epoch)))
    }
    
    data_is = data_extra_full$data
    extra_is = data_extra_full$extra
      
    # EDA checks
    if (FALSE) {
      table(extra_is$Y)
      table(extra_oos$Y)
      plot(extra_is$Y[data_is$seq_ID]+runif(data_is$R, -0.2, 0.2), data_is$Z+runif(data_is$R, -0.2, 0.2),
           col=colorspace::rainbow_hcl(data_is$A, alpha=0.8)[data_is$a], pch=16,
           xlab='Y', ylab='Z')
      legend('bottomright', legend=1:data_is$A, col=colorspace::rainbow_hcl(data_is$A),
             pch=16, horiz=T)
      plot(extra_is$Y[data_is$seq_ID]+runif(data_is$R, -0.2, 0.2), apply(data_is$C, 1, which.max)+runif(data_is$R, -0.2, 0.2),
           col=rgb(0,0,0,alpha=0.6), pch=16,
           xlab='Y', ylab='max C')
      sapply(1:data_is$L, 
             \(l) (data_is$C[(extra_is$Y==l)[data_is$seq_ID],] |> 
                     apply(1, which.max) |> table()) / 
               sum(data_is$C[(extra_is$Y==l)[data_is$seq_ID],])) |> 
        round(2)
      boxplot(mu_fn(data_is$X, extra_is$beta0, extra_is$beta) ~ ordered(extra_is$Y), horizontal=T,
              xlab='beta0 + X %*% beta', ylab='Y')
    }
    
    
    ## simulation hyperparameters ##
    
    # nice starting params using linear fit
    linfit_beta = rep(0, data_is$p)  
    
    # save full datasets
    if (!real_data) {
      save(data_extra_full, file=here::here(save_path, sprintf('data_extra_full_%i.RData', epoch)))
    }
  }
  
  ## run simulations ##
  
  # iterate through each model setting
  for (setting in 1:n_sim_settings) {
    
    # iterate through each version of the model setting
    for (version in 1:n_sim_versions[setting]) {
      
      print(paste('Running ', sim_setting_names[[setting]][version], '...', sep=''))
      
      results = NA
      
      # convert data to fit model setting #
      
      if (real_data) {
        data_extra_conv = convert_sim_data(data_extra_full, setting)
      } else {
        max_Z_labelled = max_Z_C_labelled_lst[[prop_Z_labelled_linear]]$Z
        if (setting %in% c(3,5))
          max_Z_labelled = max_Z_C_labelled_lst[[version]]$Z
        
        max_C_labelled = max_Z_C_labelled_lst[[version]]$C
        
        cur_C_threshold = NA
        if (setting %in% c(1,2)) {
          cur_C_threshold = C_threshold_lst[version]
        }
        
        data_extra_conv = convert_sim_data(data_extra_full, setting,
                                           max_Z_labelled, max_C_labelled,
                                           cur_C_threshold)
      }
      
      data = data_extra_conv$data
      extra = data_extra_conv$extra
      fallible = data_extra_conv$fallible
      
      ## MCMC samples ##
      
      if (load) {
        load(here::here(save_path, sprintf('results_%i_%i_%i.RData', setting, version, epoch)))
      } else {
      
        if (setting==1) {
          
          start_time = proc.time()
          
          results = tryCatch(
            {linear_model(data, n_iter=n_iter, n_burn=n_burn, seed=seed_epoch+1e3)},
            error={\(condition) {
              message(sprintf('Error in setting %i version %i epoch %i.', setting, version, epoch))
              message(condition)
              NA
              }})
          
          end_time = proc.time()
          
          runtime_results = c(setting, version, epoch, (end_time - start_time)[3])
          
          if (FALSE) {
            tracedensplots(results$samples, data, extra)
          }
            
          linfit_beta = colMeans(results$samples$beta)
          
        } else {
          
          label = paste(setting, version, epoch, sep='_')
          
          start_time = proc.time()
            
          # burnin
          burnin = MH_alg(n_iter=n_burn, data,
                          params0=list(beta=linfit_beta),
                          adaptive_sampling=T,
                          fallible=fallible,
                          label=label,
                          seed=seed_epoch+1e3)
  
          # MCMC results
          results = MH_alg(n_iter=n_iter, data, 
                           params0=burnin$params, 
                           sampling_params=burnin$sampling_params,
                           prior0=burnin$hp$prior,
                           prop0=burnin$hp$prop,
                           fallible=burnin$fallible,
                           adaptive_sampling=F,
                           label=label,
                           seed=seed_epoch+2e3)
  
          end_time = proc.time()
          
          runtime_results = c(setting, version, epoch, (end_time - start_time)[3])
          
          if (FALSE) {
            print_acceptance_rates(burnin$cluster_acceptance_rates)
            print_acceptance_rates(results$acceptance_rates)
            sum(burnin$params$Y!=extra$Y)
            which(burnin$params$Y!=extra$Y)
            tracedensplots(burnin$samples, data, extra, burnin=200)
            tracedensplots(burnin$samples, data, burnin=0)
            tracedensplots(results$samples, data, extra, burnin=0)
            tracedensplots(results$samples, data, burnin=0)
          }
        }
      }
      
      ## analysis ##
      
      if (calc_RPS) {
        
        Y_rps_oos = NA; Y_rps_is = NA
        if (setting!=1) {
          # Y out-of-sample
          Y_rps_oos = Y_rps_oos_fn(data_oos, extra_oos, results)

          # Y in-sample
          Y_rps_is = Y_rps_is_fn(data, extra, results)
        }
          
        # save RPS results
        RPS_results = c(setting, version, epoch, Y_rps_oos, Y_rps_is)
        
        write.table(matrix(RPS_results, nrow=1), file=here::here(save_path, 'RPS_results.csv'), sep=',', append=T, row.names=F, col.names=F)
      }
        
      if (!load) {
        # save results
        save(results, file=here::here(save_path, sprintf('results_%i_%i_%i.RData', setting, version, epoch)))
        
        write.table(matrix(runtime_results, nrow=1), file=here::here(save_path, 'runtime_results.csv'), sep=',', append=T, row.names=F, col.names=F)
      }
    }
  }

}
  
close(pb)
stopCluster(cluster) 

