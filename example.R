
################################################################################

# Cohen et al. 2026                                                 
# Example file

################################################################################


# source files
source('models/OrdinalCompModelSetup.R')
source('models/OrdinalCompSimFuncs.R')

# load first simulated dataset
load('sim_results/data_extra_full_1.RData')

# randomly select 10% of images of annotate
Z_labelled = replace(rep(F, data_extra_full$data$R), sample(1:data_extra_full$data$R, 0.1*data_extra_full$data$R, replace=F), T)

# apply annotation
data_extra_conv = convert_sim_data(data_extra_full, setting=5, Z_labelled)

# fit model
results = MH_alg(n_iter=2000, data_extra_conv$data, seed=100)

# plot trace and density plots for beta
tracedensplots(results$samples, data_extra_conv$data, data_extra_conv$extra, params='beta', plots_across=2)
