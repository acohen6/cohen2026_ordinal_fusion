
################################################################################

# Cohen et al. 2026                                                 
# Settings for simulation study

################################################################################


## imports ##

library(dplyr)


## load in model functions ##
sim_func_path = 'models/OrdinalCompSimFuncs.R'
source(sim_func_path)

## global sim settings ##

# BEGIN - NEEDED FOR MAKING FIGURES #

# what to run
load = F  # reload old data? 
calc_RPS = T  # calculate RPS?
real_data = F  # real data (or sim)?
# WARNING -- will rewrite files

# folder where sim saves are stored
if (real_data) {
  save_path = 'real_data'  
} else {
  save_path = 'sim_results'
}

n_epochs = 100  # number of replicates at each setting

C_threshold_lst = sort(c(0.75, 0.9, 0.99))  # threshold at which to filter confidences
prop_Z_labelled_lst = sort(c(0.1, 0.2, 0.5))  # proportion of images which are manually annotated
prop_Z_labelled_linear = 2  # which proportion to use for the linear fit
sim_setting_names = list(linear=paste('linear (', round(100*C_threshold_lst), '%)', sep=''), 
                         threshold=paste('maximum (', round(100*C_threshold_lst), '%)', sep=''), 
                         ordinal=paste('ordinal-only (', 100*prop_Z_labelled_lst, '%)', sep=''), 
                         compositional='compositional-only', 
                         full=paste('full model (', 100*prop_Z_labelled_lst, '%)', sep=''))
n_sim_settings = length(sim_setting_names)
n_sim_versions = sapply(sim_setting_names, length)

# END - NEEDED FOR MAKING FIGURES #

# seed for run
starting_seed = 32469

n_burn = 10000  # burn-in iterations
n_iter = 5000  # saved iterations

N_is = 500  # number of in-sample seqs
N_oos = 1000  # number of out-of-sample seqs

# covariates to sample from
covariates_all = NA
if (!load) {
  covariates_all = read.csv(here::here(save_path, 'image_data_cov_CC.csv'))
  covariates_all = covariates_all %>%
    select(is_buck_rut, wtdeer_freq, pop_sqmi, 
           landscape_diversity, green_open, green_forest,
           is_day) %>%
    na.omit()
}

