
################################################################################

# Cohen et al. 2026                                                 
# Code to analyze results of simulation study

################################################################################


## imports ##

library(latex2exp)

## setup ##

source('models/OrdinalCompModelSetup.R')

## reload saved files ##

data_is_lst = list()
extra_full = NA
converted_data_extra_lst = list()
for (epoch in 1:n_epochs) {
  load(here::here(save_path, sprintf('data_extra_full_%i.RData', epoch)))
  data_is_lst[[epoch]] = data_extra_full$data
  extra_full = data_extra_full$extra
  
  load(here::here(save_path, sprintf('max_Z_C_labelled_%i.RData', epoch)))
  converted_data_extra_lst[[epoch]] = list()
  for (setting in 1:n_sim_settings) {
    converted_data_extra_lst[[epoch]][[setting]] = list()
    for (version in 1:n_sim_versions[setting]) {
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
      
      converted_data_extra_lst[[epoch]][[setting]][[version]] = data_extra_conv
    }
  }
}

results_lst = list()
for (setting in 1:n_sim_settings) {
  results_lst[[setting]] = list()
  for (version in 1:n_sim_versions[setting]) {
    results_lst[[setting]][[version]] = list()
    for (epoch in 1:n_epochs) {
      results = NA
      tryCatch({
        suppressWarnings(load(here::here(save_path, sprintf('results_%i_%i_%i.RData', setting, version, epoch))))
      },
      error={\(condition) {
        message(sprintf('Warning: cannot load setting %i version %i epoch %i.', setting, version, epoch))
      }})
      
      results_lst[[setting]][[version]][[epoch]] = results
    }
  }
}


# pretty table function
format_sim_table = function(data, colnames, cornername='', decimals=2, bold_best=F) {
  p = ncol(data)
  
  sim_setting_names_table = list(linear=paste('\\quad Thresholded at $', round(100*C_threshold_lst), '\\%$', sep=''), 
                                 threshold=paste('\\quad Thresholded at $', round(100*C_threshold_lst), '\\%$', sep=''), 
                                 ordinal=paste('\\quad $', 100*prop_Z_labelled_lst, '\\%$ annotated', sep=''), 
                                 compositional='(4) Compositional-only', 
                                 full=paste('\\quad $', 100*prop_Z_labelled_lst, '\\%$ annotated', sep=''))
  
  sim_setting_names_table = sim_setting_names_table |> unlist()
  
  
  header = paste(paste(c(cornername, colnames), collapse=' & '), ' \\\\', collapse='')
  sections = paste(c('(1) Linear', '(2) Maximum', '(3) Ordinal-only', NA, '(5) Full model'), paste(rep('&', p), collapse=' '))
  hline = '\\hline'
  
  # format numbers
  if (!bold_best) {
    data = apply(data, 1, \(x) sprintf(paste('$%.', decimals, 'f$', sep=''), x)) |> t()
  } else {
    data = apply(data, 1, \(x) sprintf(paste('%.', decimals, 'f', sep=''), x)) |> t()
    for (col in 1:ncol(data)) {
      best = which.min(data[,col])
      data[best, col] = paste('\\mathbf{', data[best,col], '}', sep='')
    }
    data = paste('$', data, '$', sep='') |> matrix(ncol=p)
  }
  
  
  # version names
  data = cbind(sim_setting_names_table, data) |> 
    apply(1, \(x) paste(x, collapse=' & '))
  
  numbers = 
    c(sections[1], data[1:3],
      sections[2], data[4:6],
      sections[3], data[7:9],
      data[10],
      sections[5], data[11:13]) |> 
    paste(collapse=' \\\\ ')
  
  tab = 
    paste(paste('\\begin{tabular}{',  paste(c('l', rep('c', p)), collapse=' '), '}', sep=''),
          hline, header, hline, numbers, '\\\\', hline, sep=' ',
          '\\end{tabular}')
  
  cat(tab)
}



# RPS results #

RPS_results_raw = read.csv(here::here(save_path, 'RPS_results.csv'))

RPS_Y0_adj = matrix(NA, nrow=0, ncol=5)
colnames(RPS_Y0_adj) = c('setting', 'version', 'version5', 'epoch', 'RPS')

for (setting in c(2,3,4)) {
  for (epoch in 1:n_epochs) {
    for (version in 1:n_sim_versions[setting]) {
      if (setting==3) {
        RPS_Y0_5 = RPS_results_raw[RPS_results_raw$setting==5 & RPS_results_raw$version==version & RPS_results_raw$epoch==epoch,]$Y_0
        RPS_Y0_cur = RPS_results_raw[RPS_results_raw$setting==setting & RPS_results_raw$version==version & RPS_results_raw$epoch==epoch,]$Y_0
        RPS_Y0_adj = rbind(RPS_Y0_adj, c(setting, version, version, epoch, RPS_Y0_cur - RPS_Y0_5))
      } else {
        for (version5 in 1:n_sim_versions[5]) {
          RPS_Y0_5 = RPS_results_raw[RPS_results_raw$setting==5 & RPS_results_raw$version==version5 & RPS_results_raw$epoch==epoch,]$Y_0
          RPS_Y0_cur = RPS_results_raw[RPS_results_raw$setting==setting & RPS_results_raw$version==version & RPS_results_raw$epoch==epoch,]$Y_0
          RPS_Y0_adj = rbind(RPS_Y0_adj, c(setting, version, version5, epoch, RPS_Y0_cur - RPS_Y0_5))
        }
      }
    }
  }
}


# for color figures
base_colors = c('goldenrod1', 'olivedrab1', 'turquoise2', 'mediumpurple3', 'maroon4')

cols_orig = sapply(1:n_sim_settings, \(setting) rep(base_colors[setting], n_sim_versions[setting]))

first_upper = \(x) {substr(x,1,1) = toupper(substr(x,1,1)); x}


# RPS plots

dev.new(height=6, width=9, units='in', noRStudioGD=T)
par(mar=c(5,11,1,3)+0.1, cex.lab=1.2, cex.axis=1.2)
layout(matrix(c(1,1,1,2,2,2,3,3,4,4,5,5), 6, 2))

# basic Y and Y_0 boxplots
for (t in 1:2) {
  name = c('Y', 'Y_0')[t]
  pretty_name = c(TeX(paste('In-sample RPS', sep='')),
                  TeX(paste('Out-of-sample RPS', sep='')))[t]
  cur_RPS = RPS_results_raw %>% as.data.frame() %>%
    filter(setting > 1) %>% rowwise() %>%
    mutate(name=sim_setting_names[[setting]][version] %>% first_upper())
  cur_RPS$name = ordered(cur_RPS$name, levels=rev(unique(cur_RPS$name)))
  cur_RPS = cbind(cur_RPS[name], 
                  cur_RPS %>% select(-Y_0, -Y))
  cur_RPS = cur_RPS[nrow(cur_RPS):1,]
  colnames(cur_RPS)[1] = 'RPS'
  
  # get mean for each setting/version
  #cur_RPS %>% group_by(name) %>% summarize(RPS = round(median(RPS), 3))
  
  boxplot(RPS~name, data=cur_RPS, 
          col=rev(unlist(cols_orig[-1])),
          ylab='', xlab=pretty_name, las=1, horizontal=T, outline=F)
}

# adjusted Y_0 boxplots
for (cur_version5 in 1:n_sim_versions[5]) {
  cur_RPS = RPS_Y0_adj %>% as.data.frame() %>% 
    filter(version5==cur_version5) %>% rowwise() %>%
    mutate(name=sim_setting_names[[setting]][version] %>% first_upper())
  cur_RPS$name = ordered(cur_RPS$name, levels=rev(unique(cur_RPS$name)))
  cur_RPS = cur_RPS[nrow(cur_RPS):1,]
  boxplot(RPS~name, data=cur_RPS,
          col=apply(cur_RPS %>% select(setting, version) %>% distinct(),
                    1, \(row) cols_orig[[row[1]]][row[2]]),
          ylab='', xlab=TeX(paste('Relative out-of-sample RPS', sep='')),
          las=1, horizontal=T, outline=F)
  abline(v=0, lty=2, lwd=2, col=rgb(0, 0, 0, alpha=0.5))
  legend('bottomright', bty='n', cex=1.2,
         legend=paste('Baseline: ', sim_setting_names[[5]][cur_version5], sep=''))
}


par(mfrow=c(1,1), mar=c(5,5,3,1)+0.1)


# MSE #

MSE = list()

for (setting in 1:n_sim_settings) {
  MSE[[setting]] = list()
  for (version in 1:n_sim_versions[setting]) {
    MSE[[setting]][[version]] = rep(NA, data_is_lst[[1]]$p)
    for (f in 1:data_is_lst[[1]]$p) {
      beta_mns = sapply(results_lst[[setting]][[version]],
                        \(epoch) mean(epoch$samples$beta[,f]))
      MSE[[setting]][[version]][f] = 
        mean((beta_mns - extra_full$beta[f])^2)
    }
  }
}

MSE_table = MSE |> unlist() |> matrix(ncol=data_is_lst[[1]]$p, byrow=T)


# coverage #

coverage = list()  # percent of CIs that contain the true value
detection = list()  # percent of CIs that contain zero

for (setting in 1:n_sim_settings) {
  coverage[[setting]] = list()
  detection[[setting]] = list()
  for (version in 1:n_sim_versions[setting]) {
    coverage[[setting]][[version]] = rep(NA, data_is_lst[[1]]$p)
    detection[[setting]][[version]] = rep(NA, data_is_lst[[1]]$p)
    for (f in 1:data_is_lst[[1]]$p) {
      coverage[[setting]][[version]][f] = 
        sapply(1:n_epochs,
               \(epoch) between(extra_full$beta[f],
                                quantile(results_lst[[setting]][[version]][[epoch]]$samples$beta[,f], 0.025),
                                quantile(results_lst[[setting]][[version]][[epoch]]$samples$beta[,f], 0.975))) |> 
        mean()
      
      detection[[setting]][[version]][f] = 
        sapply(1:n_epochs,
               \(epoch) between(0,
                                quantile(results_lst[[setting]][[version]][[epoch]]$samples$beta[,f], 0.025),
                                quantile(results_lst[[setting]][[version]][[epoch]]$samples$beta[,f], 0.975))) |> 
        mean()
    }
  }
}

coverage_table = coverage |> unlist() |> matrix(ncol=data_is_lst[[1]]$p, byrow=T)

detection_table = detection |> unlist() |> matrix(ncol=data_is_lst[[1]]$p, byrow=T)


# latex format for MSE, coverage, detection

beta_equals_names = sapply(1:data_is_lst[[1]]$p, \(x) sprintf('$\\beta_{%i}=%.1f$', x, extra_full$beta[x]))

format_sim_table(1e3*MSE_table, beta_equals_names, cornername='Scenario', decimals=2, bold_best=T)

format_sim_table(coverage_table, beta_equals_names, cornername='Scenario', decimals=2)

format_sim_table(cbind(1-detection_table[,extra_full$beta!=0],
                       detection_table[,extra_full$beta==0]), 
                 c(beta_equals_names[extra_full$beta!=0], beta_equals_names[extra_full$beta==0]),
                 cornername='Scenario',
                 decimals=2)



# quantiles #

CIs = array(NA, dim=c(n_sim_settings, max(n_sim_versions), n_epochs, data_is_lst[[1]]$p, 3), 
            dimnames=list(setting=sim_setting_names,
                          version=1:max(n_sim_versions),
                          epoch=1:n_epochs,
                          beta=paste('beta', 1:data_is_lst[[1]]$p, sep='_'), 
                          quantile=paste('quant', c(2.5, 50, 97.5), sep='_')))

for (setting in 1:n_sim_settings) {
  for (version in 1:n_sim_versions[setting]) {
    for (epoch in 1:n_epochs) {
      for (f in 1:data_is_lst[[1]]$p) {
        CIs[setting,version,epoch,f,] = quantile(results_lst[[setting]][[version]][[epoch]]$samples$beta[,f], c(0.025, 0.5, 0.975))
      }
    }
  }
}


# valid #

valid = data.frame(matrix(nrow=0, ncol=5))
for (epoch in 1:n_epochs) {
  for (setting in 1:n_sim_settings) {
    for (version in 1:n_sim_versions[setting]) {
      valid = rbind(valid, 
                    c(epoch, setting, version,
                      converted_data_extra_lst[[epoch]][[setting]][[version]]$data$img_used |> sum(),
                      converted_data_extra_lst[[epoch]][[setting]][[version]]$data$N))
    }
  }
}
colnames(valid) = c('epoch', 'setting', 'version', 'images', 'sequences')

num_valid = valid %>% group_by(setting, version) %>%
  select(images, sequences) %>%
  summarize(across(.cols=c(images, sequences), median)) %>% 
  ungroup() %>%
  arrange(setting, version)

format_sim_table(data=num_valid %>% select(images, sequences), colnames=c('Images', 'Sequences'), cornername='Scenario', decimals=0)

# runtime results #

runtime_results = read.csv(paste(save_path, 'runtime_results.csv', sep=''))
runtime_results = runtime_results %>% group_by(setting) %>% 
  summarize(across(.cols=-1, .fns=\(x) mean(x)/60))

runtime_results

knitr::kable(runtime_results, format='latex', digits=0)


# plot CI's with true value

cols = sapply(cols_orig, colorspace::adjust_transparency, alpha=1/n_epochs)
within_sp = 0.1
btwn_sp = 0.4
within_ticks = seq(0, 4*within_sp, by=within_sp)
cex = 0.8

dev.new(height=4, width=9, units='in', noRStudioGD=T)
par(mfrow=c(1,1), mar=c(3,4,1,10), oma=c(0,0,0,0))
plot(0, xlim=c(0, (btwn_sp+4*within_sp)*data_is_lst[[1]]$p), ylim=range(CIs, na.rm=T), 
     main='', xlab='', ylab='posterior',
     xaxt='n', yaxt='n', bty='l', type='n', cex.lab=cex)
axis(1, at=seq(2*within_sp, (btwn_sp+4*within_sp)*data_is_lst[[1]]$p, by=btwn_sp+4*within_sp),
     #labels=1:data_is$p, 
     labels=sapply(1:data_is_lst[[1]]$p, \(x) TeX(sprintf(r'($\beta_{%i}$)', x))),
     cex.axis=cex)
axis(2, cex.axis=cex)
legend('right', legend=sim_setting_names, fill=cols_orig, cex=cex,
       inset=c(-0.35, 0), xpd=T, bty='n')
abline(h=0, lty=2, lwd=2)
for (f in 1:data_is_lst[[1]]$p) {
  start_tick = (btwn_sp+4*within_sp)*(f-1)
  arrows(x0=start_tick-within_sp, x1=start_tick+4*within_sp+within_sp,
         y0=extra_full$beta[f], y1=extra_full$beta[f],
         length=0, lwd=2)
  for (epoch in 1:n_epochs) {
    arrows(x0=start_tick+within_ticks, x1=start_tick+within_ticks, 
           y0=CIs[,epoch,f,1], y1=CIs[,epoch,f,3], 
           angle=90, lwd=4, length=0, col=cols)
  }
}

