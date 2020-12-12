library(ggplot2)
library(reshape2)
library(plotly)

setwd("~/Documents/Research/Climate/Data")
all_files <- list.files(pattern = '*.dat', path = 'CMIP_data')
string_part1 <- sapply(all_files, function(x) {return (strsplit(x, '_')[[1]][1])})
string_part2 <- sapply(all_files, function(x) {return (strsplit(x, '_')[[1]][2])})
string_part3 <- sapply(all_files, function(x) {return (strsplit(x, '_')[[1]][3])})
string_part4 <- sapply(all_files, function(x) {return (strsplit(x, '_')[[1]][4])})
string_part5 <- sapply(all_files, function(x) {return (strsplit(x, '_')[[1]][5])})
string_part6 <- sapply(all_files, function(x) {return (strsplit(x, '_')[[1]][6])})
all_files_split <- data.frame(v1 = string_part1, v2 = string_part2, v3 = string_part3,
                              v4 = string_part4, v5 = string_part5, v6 = string_part6)

# temp <- read.table('global_tas_Amon_bcc-csm1-1_historicalNat_r1i1p1.dat', header = F)
# temp1 <- read.table('global_tas_Amon_bcc-csm1-1_rcp26_r1i1p1.dat', header = F)
# temp2 <- read.table('global_tas_Amon_bcc-csm1-1_rcp45_r1i1p1.dat', header = F)

temp <- read.table('had4_krig_annual_v2_0_0.txt')

all_data <- list()
count <- 0
for (file in all_files)
{
  temp <- read.table(paste0('CMIP_Data/',file), header = F)
  count <- count + 1
  if (all_files_split$v5[count] == 'historicalNat') {next}
  all_data[[paste0('_s_', all_files_split$v4[count], '_s_', all_files_split$v5[count], 
                   '_s_', all_files_split$v6[count])]] <- temp
  print(dim(temp))
  
}

max_year <- 0
min_year <- 2000
for (data_set in all_data)
{
  max_year <- max(max_year, data_set$V1)
  min_year <- min(min_year, data_set$V1)
}

cleaned_all_data <- list()
months_average <- as.character(months(as.Date('2017-01-01') + 31 * 0:11))
months_average[13] <- 'Average'
for (i in 1:12)
{
  cleaned_all_data[[months_average[i]]] <- data.frame(year = min_year:max_year)
  count <- 2
  for (model in names(all_data))
  {
    cleaned_all_data[[months_average[i]]] <- merge(cleaned_all_data[[months_average[i]]], 
                                                   all_data[[model]][, c('V1', paste0('V', i+1))], 
                                                   by.x = 'year', 
                                                   by.y = 'V1', all.x = TRUE)
    names(cleaned_all_data[[months_average[i]]])[count] <- model
    if (any(cleaned_all_data[[months_average[i]]][, count] < 273, na.rm = T))
    {
      cleaned_all_data[[months_average[i]]][which(cleaned_all_data[[months_average[i]]][, count] < 273), count] <- NA
    }
    count <- count + 1
  }
  cleaned_all_data[[months_average[i]]]$Quant025 <- apply(cleaned_all_data[[months_average[i]]][,-1], 1, quantile, probs = 0.025, na.rm = TRUE)
  cleaned_all_data[[months_average[i]]]$Quant050 <- apply(cleaned_all_data[[months_average[i]]][,-1], 1, quantile, probs = 0.05, na.rm = TRUE)
  cleaned_all_data[[months_average[i]]]$Quant500 <- apply(cleaned_all_data[[months_average[i]]][,-1], 1, quantile, probs = 0.5, na.rm = TRUE)
  cleaned_all_data[[months_average[i]]]$Quant950 <- apply(cleaned_all_data[[months_average[i]]][,-1], 1, quantile, probs = 0.95, na.rm = TRUE)
  cleaned_all_data[[months_average[i]]]$Quant975 <- apply(cleaned_all_data[[months_average[i]]][,-1], 1, quantile, probs = 0.975, na.rm = TRUE)
  cat(months_average[i], '\n')
}

cleaned_all_data$Average <- cleaned_all_data[[1]][, 1:344]
for (i in 2:12)
{
  cleaned_all_data$Average <- cleaned_all_data$Average + cleaned_all_data[[i]][, 1:344]
}
cleaned_all_data$Average <- cleaned_all_data$Average / 12
cleaned_all_data$Average$Quant025 <- apply(cleaned_all_data$Average[,-1], 1, quantile, probs = 0.025, na.rm = TRUE)
cleaned_all_data$Average$Quant050 <- apply(cleaned_all_data$Average[,-1], 1, quantile, probs = 0.05, na.rm = TRUE)
cleaned_all_data$Average$Quant500 <- apply(cleaned_all_data$Average[,-1], 1, quantile, probs = 0.5, na.rm = TRUE)
cleaned_all_data$Average$Quant950 <- apply(cleaned_all_data$Average[,-1], 1, quantile, probs = 0.95, na.rm = TRUE)
cleaned_all_data$Average$Quant975 <- apply(cleaned_all_data$Average[,-1], 1, quantile, probs = 0.975, na.rm = TRUE)

cleaned_all_data$Average_Diff <- cleaned_all_data$Average[, 1:344]
all_models_split <- unlist(strsplit(names(cleaned_all_data$Average_Diff), '_s_'))[-1]
all_models <- unique(all_models_split[seq(2, 1370, 4)])
selected_indices <- c()
for (model in all_models)
{
  all_trajectories <- all_models_split[which(all_models_split == model) + 2]
  all_scenarios <- all_models_split[which(all_models_split == model) + 1]
  indices <- (which(all_models_split == model) + 2) / 4
  if ('ave.dat' %in% all_trajectories)
  {
    if (!is.na(which(all_trajectories == 'ave.dat')[1]))
    {
      selected_indices <- c(selected_indices, indices[which(all_trajectories == 'ave.dat')])
    }
    else
    {
      selected_indices <- c(selected_indices, indices[which(all_trajectories == 'ave.dat')])
    }
  }
  else
  {
    selected_indices <- c(selected_indices, indices)
  }
}

all_df <- data.frame()
for (index in (selected_indices + 1))
{
  temp_df <- cleaned_all_data$Average_Diff[, c(1, index)]
  names(temp_df) <- c('year', 'temperature')
  pre_industrial_avg <- mean(subset(temp_df, year >= 1981 & year <= 2005)$temperature)
  temp_df$temperature <- temp_df$temperature - pre_industrial_avg
  temp_df <- subset(temp_df, year >= 2006 & year <= 2100)
  val <- strsplit(names(cleaned_all_data$Average_Diff)[index], split = '_s_')[[1]]
  temp_df$model <- val[2]
  temp_df$scenario <- val[3]
  temp_df$co2 <- rcp.carbon.cum.complete.adjusted[, val[3]]
  all_df <- rbind(all_df, temp_df)
}

r2 <- c()
slope <- c()
for (model_x in unique(all_df$model))
{
  temp_df <- subset(all_df, model == model_x)
  m1 <- lm(temperature ~ co2, data = temp_df)
  r2 <- c(r2, summary(m1)$r.squared)
  slope <- c(slope, summary(m1)$coef[2,1])
}

all_hist_df <- data.frame()

model <- 'initial'
for (index in (selected_indices + 1))
{
  temp_df <- cleaned_all_data$Average_Diff[, c(1, index)]
  val <- strsplit(names(cleaned_all_data$Average_Diff)[index], split = '_s_')[[1]]
  model_x <- val[2]
  if (model_x == model)
  {
    next
  }
  else
  {
    model <- model_x
  }
  names(temp_df) <- c('year', 'temperature')
  pre_industrial_avg <- mean(subset(temp_df, year >= 1981 & year <= 2005)$temperature)
  temp_df$temperature <- temp_df$temperature - pre_industrial_avg
  temp_df <- subset(temp_df, year >= 1861 & year <= 2005)
  temp_df$model <- val[2]
  temp_df$real_data <- real_data[12:156, 'V6'] - mean(real_data[132:156, 'V6'])
  all_hist_df <- rbind(all_hist_df, temp_df)
}

for (i in 1:39)
{
  pdf(file = paste0('CMIP5_models/', all_models[i], '.pdf'), width = 8.5, height = 5)
  
  p <- ggplot(subset(all_hist_df, model == all_models[i]), aes(x = year)) + 
    geom_line(aes(y=temperature), alpha=1, size=1.3, color="#e41a1c") + 
    geom_line(aes(y=real_data), alpha=1, size=1.3, color="black") + 
    ggtitle(all_models[i]) + ylab('Anomaly')
  gridExtra::grid.arrange(p)
  dev.off()
}

make_bugs_data <- function(all_hist_df, model_x)
{
  year <- 1861:2005
  temp_df <- subset(all_hist_df, model==model_x)
  model_forecast <- temp_df$temperature
  hadCrut_observation <- temp_df$real_data
  sd_hadCrut <- real_data$V3[real_data$V1 %in% year]
  model_forecast_anomaly <- model_forecast
  observed_anomaly <- hadCrut_observation
  n_years <- length(year)
  
  return (list(n_years = n_years, 
               sd_delta = sd_hadCrut,
               model_forecast = model_forecast_anomaly,
               observed_anomaly = observed_anomaly))
}

fit_ar1_model <- function(all_hist_df, model_x,
                          n.iterations = 11000, n.adapt = 1000, n.chains = 3, thin = 1, 
                          model_name = 'bayes_ar1.bug',
                          var.list = c("true_anomaly", "rho", "sd_w"))
{
  library(rjags)
  library(coda)
  library(reshape2)
  print("Calling make.bugs.data")
  jags.input <- make_bugs_data(all_hist_df, model_x)
  # browser()
  print("Successfully called make.bugs.data")
  jags.out <- jags.model(model_name, 
                         data=jags.input,
                         n.chains=n.chains, n.adapt=n.adapt)
  jags.output <- coda.samples(jags.out,
                              var.list,
                              n.iter=n.iterations,
                              thin=thin)
  return (list(input = jags.input, output = jags.output))
  
}

fit_result <- fit_ar1_model(all_hist_df, "MRI-CGCM3", n.iterations = 10000)

true_anomaly_names <- paste('true_anomaly[', 1:145, ']', sep = '')
posterior_anomaly <- as.matrix(fit_result[[2]][[1]][, true_anomaly_names])
summary_anomaly <- data.frame(year = 1861:2005)
summary_anomaly$median <- apply(posterior_anomaly, 2, median)
summary_anomaly$q025 <- apply(posterior_anomaly, 2, quantile, probs = 0.025)
summary_anomaly$q100 <- apply(posterior_anomaly, 2, quantile, probs = 0.1)
summary_anomaly$q900 <- apply(posterior_anomaly, 2, quantile, probs = 0.9)
summary_anomaly$q975 <- apply(posterior_anomaly, 2, quantile, probs = 0.975)
summary_anomaly$observed_anomaly <- fit_result$input$observed_anomaly

ggplot(data = summary_anomaly, aes(x = year))  +
  geom_line(aes(y = median, color = 'True'), size = 2) + scale_color_manual(values = c('HadCrut4' = 'black', 'True' = 'red')) + 
  geom_line(aes(y = observed_anomaly, color = 'HadCrut4')) + 
  geom_ribbon(aes(ymin = q025, ymax = q975), fill = 'red', alpha = 0.2) + 
  geom_ribbon(aes(ymin = q100, ymax = q900), fill = 'red', alpha = 0.3) + xlab('year') + ylab('Anomaly')


project_xt <- function(carbon_trajectory, rcp_co2, model_x, year.range = c(2005, 2110), start.from.2015=TRUE)
{
  co2.trajectory <- carbon_trajectory$CO2 / 1e9
  carbon.cum <- 0
  proj_length <- length(year.range[1]:year.range[2]) - 1
  co2.trajectory.combined <- co2.trajectory[-1]
  if (start.from.2015)
  {
    for (year in 2015:2006)
    {
      co2.trajectory.combined <- c(sum(data.medium$CO2.total[data.medium$Year == year] * 11/3)/1000, co2.trajectory.combined)
    }
  }
  cum_co2 <- cumsum(co2.trajectory.combined)
  # cum_co2 <- co2.trajectory[1] * 5 + cum_co2
  rcp_co2_adjusted <- rcp_co2
  rcp_co2_adjusted <- rcp_co2_adjusted[c(3,2,1,4), ]
  rcp_co2_adjusted <- t(as.matrix(rcp_co2_adjusted[, -(1:2)]))
  rcp_co2_adjusted <- rcp_co2_adjusted[1:(length(year.range[1]:year.range[2]) - 1), ]
  rcp_co2_adjusted <- as.data.frame(rcp_co2_adjusted)
  names(rcp_co2_adjusted) <- paste0('rcp', c(26, 45, 60, 85))
  rcp_co2_adjusted$year <- (year.range[1]+1):year.range[2]
  rcp_temp_adjusted <- subset(all_df, model == model_x)
  temp1 <- rcp_temp_adjusted$temperature[1]
  rcp_temp_adjusted$temperature <- rcp_temp_adjusted$temperature - temp1
  temp_model <- lm(temperature ~ co2 + 0, data = rcp_temp_adjusted)
  
  # projected_temp <- mean(as.numeric(rcp_temp_adjusted[1, -1]))
  projected_temp <- as.numeric(predict(temp_model, newdata = data.frame(co2=cum_co2))) + temp1
  return (projected_temp)
}

project_zt <- function(fit_result, index, year.end = 2100, year.start = 2005, proj_xt)
{
  fit_result_chain <- fit_result$output[[1]]
  zt2005 <- fit_result_chain[index, 'true_anomaly[145]']
  rho <- fit_result_chain[index, 'rho']
  sd_w <- fit_result_chain[index, 'sd_w']
  # D <- fit_result_chain[index, 'D']
  D <- 0
  zt_proj <- numeric(year.end - year.start + 1)
  if (year.start == 2005)
  {
    zt2005 <- zt2005 - fit_result$input$model_forecast[145]
    zt_proj[1] <- zt2005
  }
  else
  {
    idx <- which(real_data_adjusted$V1 == year.start)
    zt_start <- rnorm(1, real_data_adjusted$V6[idx], real_data_adjusted$V3[idx])
    zt_proj[1] <- zt_start - proj_xt[year.start - 2005]
  }
  for (i in 2:(year.end - year.start + 1))
  {
    zt_proj[i] <- (zt_proj[i - 1] - D) * rho + D + rnorm(1, 0, sd_w)
  }
  return (zt_proj[-1])
}

project_temperature <- function(fit_result_chain, carbon.projection.trajectories, rcp_co2, model_x, n.traj = 1000, 
                                year.range = c(2005, 2110), ratio_co2 = 1, start.from.2015=TRUE)
{
  num_carbon_trajectories <- length(carbon.projection.trajectories)
  num_zt_trajectories <- nrow(fit_result_chain$output[[1]])
  carbon_trajs <- sample(num_carbon_trajectories, n.traj, replace = TRUE)
  zt_trajs <- sample(num_zt_trajectories, n.traj, replace = TRUE)
  length_proj <- year.range[2] - year.range[1]
  proj_temps <- matrix(nrow = length_proj, ncol = n.traj + 1)
  proj_xts <- matrix(nrow = length_proj, ncol = n.traj + 1)
  proj_zts <- matrix(nrow = length_proj, ncol = n.traj + 1)
  proj_temps[, 1] <- proj_zts[, 1] <- proj_xts[, 1] <- (year.range[1]+1):year.range[2]
  mean1861_1880 <- numeric(n.traj)
  for (i in 1:n.traj)
  {
    proj_xt <- project_xt(carbon.projection.trajectories[[carbon_trajs[i]]] * ratio_co2, rcp_co2, 
                          model_x, year.range, start.from.2015=start.from.2015)
    proj_zt <- project_zt(fit_result_chain, zt_trajs[i], year.end = year.range[2], year.start = year.range[1], proj_xt=proj_xt)
    proj_xt <- proj_xt[((year.range[1]+1):year.range[2])-2005]
    # proj_zt_cutting <- proj_zt[5 * (1:length_proj)]
    # proj_zt_cutting <- proj_zt[seq(year.range[1], year.range[2], 5) - 2005]
    proj_temp <- proj_xt + proj_zt # proj_zt_cutting
    proj_temps[, i + 1] <- proj_temp
    proj_xts[, i + 1] <- proj_xt
    proj_zts[, i + 1] <- proj_zt
  }
  
  return(list(projection = proj_temps, xt_projection = proj_xts, zt_projection = proj_zts))
}

calc_likelihood <- function(fit_result, idx)
{
  posterior_anomaly <- as.matrix(fit_result[[2]][[1]][, true_anomaly_names])
  z_it <- as.numeric(posterior_anomaly[idx, ])
  y_it <- as.numeric(real_data$V6 - mean(subset(real_data, V1 %in% 1981:2005)$V6))[12:156]
  V_t <- real_data$V3[12:156]
  sigma <- fit_result[[2]][[1]][idx, 'sd_w']
  print (sigma)
  rho <- fit_result[[2]][[1]][idx, 'rho']
  likelihood <- prod(1/sqrt(2 * pi)/V_t * exp(-1/2/V_t^2 * (y_it - z_it)^2))
  likelihood2 <- prod(1/sqrt(2 * pi)/sigma * exp(-1/2/sigma^2 * (z_it[2:145] - rho * z_it[1:144])^2))
  return (likelihood * likelihood2)
}

all_fits <- list()
all_projections <- list()
combined_projections <- c()

real_data_adjusted <- real_data
real_data_adjusted$V6 <- real_data_adjusted$V6 - mean(real_data_adjusted$V6[132:156])

for (i in 1:39)
{
  model_x <- all_models[i]
  fit_result <- fit_ar1_model(all_hist_df, model_x, n.iterations = 1000)
  all_fits[[model_x]] <- fit_result
}

proj_temp_res <- c()
probs <- c()
for (j in 1:4)
{
  if (j == 1) tmp <- proj.evals.2015.ar1.const$trajs.annual.worldwide
  else if (j == 2) tmp <- proj.evals.2015.adjusted$trajs.annual.worldwide
  else if (j == 3) tmp <- proj.evals.2015.adjusted$trajs.annual.worldwide.cont
  else tmp <- proj.evals.2015.adjusted$trajs.annual.worldwide.usa
  
  combined_projections <- c()
  for(i in 1:39)
  {
    model_x <- all_models[i]
    proj_list <- project_temperature(all_fits[[model_x]], tmp, rcp.carbon.cum.complete, 
                                     model_x, year.range = c(2015, 2100), n.traj = 100)
    cat(model_x, ' finished\n')
    all_projections[[model_x]] <- proj_list
    combined_projections <- cbind(combined_projections, all_projections[[model_x]]$projection[,2:101])
  }
  probs <- c(probs, mean(combined_projections[85,] < (1.5 + mean(real_data_adjusted$V6[12:31])-mean(real_data_adjusted$V6[132:156]))))
  combined.results <- t(apply(combined_projections, 1, quantile, probs=c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95)))
  combined.results <- as.data.frame(combined.results)
  combined.results <- combined.results + mean(real_data_adjusted$V6[132:156]) - mean(real_data_adjusted$V6[12:31])
  combined.results$Year <- 2016:2100
  proj_temp_res <- rbind(proj_temp_res, combined.results)
}

proj_temp_res$adjusted <- NA
proj_temp_res$adjusted[1:85] <- 'None'
proj_temp_res$adjusted[86:170] <- 'Adjusted'
proj_temp_res$adjusted[171:255] <- 'Continued'
proj_temp_res$adjusted[256:340] <- 'USA Excluded'
proj_temp_res <- proj_temp_res[, c("q025", "q050", "q100", "median", "q900", "q950", "q975", "year", "adjusted")]
names(proj_temp_res) <- c("q025", "q050", "q100", "median", "q900", "q950", "q975", "year", "adjusted")
summary_anomaly[, 2:9] <- summary_anomaly[, 2:9] - mean(summary_anomaly$observed_anomaly[1:20])
proj_temp_res$adjusted <- factor(proj_temp_res$adjusted, levels = c('None', 'Adjusted', 'Continued', 'USA Excluded'))
plot_temp_new(proj_temp_res, summary_anomaly)


combined_projections <- c()
for(i in 1:39)
{
  model_x <- all_models[i]
  proj_list <- project_temperature(all_fits[[model_x]], proj.evals.2015.ar1.const$trajs.annual.worldwide, ratio_co2 = 0.1, rcp.carbon.cum.complete, 
                                   model_x, year.range = c(2015, 2100), n.traj = 100)
  cat(model_x, ' finished\n')
  all_projections[[model_x]] <- proj_list
  combined_projections <- cbind(combined_projections, all_projections[[model_x]]$projection[,2:101])
}
combined.results <- t(apply(combined_projections, 1, quantile, probs=c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95)))
combined.results <- as.data.frame(combined.results)
combined.results <- combined.results + mean(real_data_adjusted$V6[132:156]) - mean(real_data_adjusted$V6[12:31])
combined.results$Year <- 2016:2100

bi_section <- function(objective = 2, prob = 0.5)
{
  obj <- 2
  ratio_co2_r <- ratio_co2 <- 1
  ratio_co2_l <- 0
  count <- 1
  while (abs(obj) > 1e-3)
  {
    ratio_co2 <- ratio_co2_l / 2 + ratio_co2_r / 2
    combined_projections <- c()
    for(i in 1:39)
    {
      model_x <- all_models[i]
      proj_list <- project_temperature(all_fits[[model_x]], proj.evals.2015.adjusted$trajs.annual.worldwide.cont, rcp.carbon.cum.complete, 
                                       model_x, year.range = c(2015, 2100), n.traj = 100, ratio_co2 = ratio_co2)
      all_projections[[model_x]] <- proj_list
      combined_projections <- cbind(combined_projections, all_projections[[model_x]]$projection[,2:101])
    }
    
    proj_temp <- proj_list$projection
    summary_proj <- data.frame(year = proj_temp[,1])
    # browser()
    summary_proj$quantile <- apply(combined_projections, 1, quantile, prob = prob)
    summary_proj$quantile <- summary_proj$quantile + mean(real_data_adjusted$V6[132:156]) - mean(real_data_adjusted$V6[12:31])
    global_warming <- summary_proj$quantile[summary_proj$year == 2100]
    obj <- global_warming - objective
    if (obj > 0)
    {
      ratio_co2_r <- ratio_co2
    }
    else
    {
      ratio_co2_l <- ratio_co2
    }
    cat('Iteration: ', count, ' Objective: ', obj, 'Ratio: ', ratio_co2, '\n')
    count <- count + 1
    if (count > 15) {break}
  }
  return (ratio_co2)
}

bi_section(objective = 2, prob=0.9)
bi_section(objective = 1.5, prob=0.9)

model.results <- t(apply(all_projections$`ACCESS1-0`$projection, 1, quantile, probs=c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95)))
model.results <- as.data.frame(model.results)
model.results$Year <- 2016:2100

names(model.results)[1:7] <- paste0('Quant', c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95))

model.results[, 1:7] <- model.results[, 1:7] + mean(real_data$V6[132:156]) - mean(real_data$V6[12:31]) 


font <- 'Times'
plot.obj <- ggplot() +
  labs(color="Scenario") +
  theme(panel.background=element_rect(fill="white"),
        panel.grid.minor=element_line(color="white"),
        panel.grid.major=element_line(color="grey90"),
        text=element_text(family=font, size=22))

color.us <- "#e41a1c"
plot.obj <- plot.obj +
  ylab('Anomaly') +
  geom_ribbon(aes(x=Year, ymin=Quant0.05, ymax=Quant0.95), data=model.results,
              alpha=0.2, fill='blue') +
  geom_ribbon(aes(x=Year, ymin=Quant0.1, ymax=Quant0.9), data=model.results,
              alpha=0.3, fill='blue') +
  geom_line(data=model.results, aes(x=Year, y=Quant0.5),
            alpha=1, size=1.3, color='blue') +
  # scale_x_continuous(breaks=c(2010, seq(2025, year.end, 25)), expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0))

real_data_adjusted <- real_data
real_data_adjusted$V6 <- real_data_adjusted$V6 - mean(real_data_adjusted$V6[12:31])

plot.obj <- plot.obj + 
  geom_line(aes(x = V1, y = V6), data = real_data_adjusted)

plot.obj <- plot.obj + ggtitle('ACCESS1-0 Mean Temperature Anomaly Forecast')
plot.obj


xt.results <- t(apply(all_projections$`ACCESS1-0`$xt_projection, 1, quantile, probs=c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95)))
xt.results <- as.data.frame(xt.results)
xt.results$Year <- 2016:2100

names(xt.results)[1:7] <- paste0('Quant', c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95))

xt.results[, 1:7] <- xt.results[, 1:7] + mean(real_data$V6[132:156]) - mean(real_data$V6[12:31]) 

zt.results <- t(apply(all_projections$`ACCESS1-0`$zt_projection, 1, quantile, probs=c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95)))
zt.results <- as.data.frame(zt.results)
zt.results$Year <- 2016:2100

names(zt.results)[1:7] <- paste0('Quant', c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95))

color.us <- "#e41a1c"
names_all <- c('x_t', 'z_t')

font <- 'Times'
plot.obj <- ggplot() +
  theme(panel.background=element_rect(fill="white"),
        panel.grid.minor=element_line(color="white"),
        panel.grid.major=element_line(color="grey90"),
        text=element_text(family=font, size=22))

plot.obj <- plot.obj +
  ylab('Anomaly') +
  geom_ribbon(aes(x=Year, ymin=Quant0.05, ymax=Quant0.95, fill='x_t'), data=xt.results,
              alpha=0.2) +
  geom_ribbon(aes(x=Year, ymin=Quant0.1, ymax=Quant0.9, fill='x_t'), data=xt.results,
              alpha=0.3) +
  geom_line(data=xt.results, aes(x=Year, y=Quant0.5, color='x_t'),
            alpha=1, size=1.3)
  # scale_x_continuous(breaks=c(2010, seq(2025, year.end, 25)), expand = c(0,0)) +
  # scale_y_continuous(expand = c(0,0))

plot.obj <- plot.obj +
  geom_ribbon(aes(x=Year, ymin=Quant0.05, ymax=Quant0.95, fill='z_t'), data=zt.results,
              alpha=0.2) +
  geom_ribbon(aes(x=Year, ymin=Quant0.1, ymax=Quant0.9, fill='z_t'), data=zt.results,
              alpha=0.3) +
  geom_line(data=zt.results, aes(x=Year, y=Quant0.5, color='z_t'),
            alpha=1, size=1.3)
  
colors <- c('blue', color.us)
names(colors) <- names_all
plot.obj <- plot.obj + ggtitle('ACCESS1-0 Mean Temperature Anomaly Forecast') + 
  scale_color_manual(values=colors) + scale_fill_manual(values=colors)
plot.obj


combined_projections <- c()
for(i in 1:39)
{
  model_x <- all_models[i]
  combined_projections <- cbind(combined_projections, all_projections[[model_x]]$projection)
}

combined.results <- t(apply(combined_projections, 1, quantile, probs=c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95)))
combined.results <- as.data.frame(combined.results)
combined.results$Year <- 2016:2100

names(combined.results)[1:7] <- paste0('Quant', c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95))

combined.results[, 1:7] <- combined.results[, 1:7] + mean(real_data$V6[132:156]) - mean(real_data$V6[12:31]) 

font <- 'Times'
plot.obj <- ggplot() +
  labs(color="Scenario") +
  theme(panel.background=element_rect(fill="white"),
        panel.grid.minor=element_line(color="white"),
        panel.grid.major=element_line(color="grey90"),
        text=element_text(family=font, size=22))

color.us <- "#e41a1c"
plot.obj <- plot.obj +
  ylab('Anomaly') +
  geom_ribbon(aes(x=Year, ymin=Quant0.05, ymax=Quant0.95), data=combined.results,
              alpha=0.2, fill=color.us) +
  geom_ribbon(aes(x=Year, ymin=Quant0.1, ymax=Quant0.9), data=combined.results,
              alpha=0.3, fill=color.us) +
  geom_line(data=combined.results, aes(x=Year, y=Quant0.5),
            alpha=1, size=1.3, color=color.us) +
  # scale_x_continuous(breaks=c(2010, seq(2025, year.end, 25)), expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0))

real_data_adjusted <- real_data
real_data_adjusted$V6 <- real_data_adjusted$V6 - mean(real_data_adjusted$V6[132:156])

plot.obj <- plot.obj + 
  geom_line(aes(x = V1, y = V6), data = real_data_adjusted)

plot.obj <- plot.obj + ggtitle('Assembled Global Mean Temperature Forecast')
plot.obj


for (i in 1:39)
{
  model_x <- all_models[i]
  fit_result <- fit_ar1_model(all_hist_df, model_x, n.iterations = 10000)
  
  true_anomaly_names <- paste('true_anomaly[', 1:145, ']', sep = '')
  posterior_anomaly <- as.matrix(fit_result[[2]][[1]][, true_anomaly_names])
  summary_anomaly <- data.frame(year = 1861:2005)
  summary_anomaly$median <- apply(posterior_anomaly, 2, median)
  summary_anomaly$q025 <- apply(posterior_anomaly, 2, quantile, probs = 0.025)
  summary_anomaly$q100 <- apply(posterior_anomaly, 2, quantile, probs = 0.1)
  summary_anomaly$q900 <- apply(posterior_anomaly, 2, quantile, probs = 0.9)
  summary_anomaly$q975 <- apply(posterior_anomaly, 2, quantile, probs = 0.975)
  summary_anomaly$observed_anomaly <- fit_result$input$observed_anomaly
  
  
  proj_list <- project_temperature(fit_result, proj.evals.2015.ar1.const$trajs.annual.worldwide, rcp.carbon.cum.complete, 
                                   model_x, year.range = c(2015, 2100), n.traj = 1000)
  proj_temp <- proj_list$projection
  proj_xts <- proj_list$xt_projection
  proj_zts <- proj_list$zt_projection
  
  summary_proj <- data.frame(year = proj_temp[,1])
  summary_proj$q025 <- apply(proj_temp[,-1], 1, quantile, probs = 0.025)
  summary_proj$q100 <- apply(proj_temp[,-1], 1, quantile, probs = 0.1)
  summary_proj$q900 <- apply(proj_temp[,-1], 1, quantile, probs = 0.9)
  summary_proj$q975 <- apply(proj_temp[,-1], 1, quantile, probs = 0.975)
  summary_proj$median <- apply(proj_temp[,-1], 1, median)
  
  p <- ggplot(data = summary_anomaly, aes(x = year))  +
    geom_line(aes(y = median, color = 'True'), size = 2) + 
    geom_line(data = summary_proj, mapping = aes(x = year, y = median, color = 'Projection')) + 
    geom_line(data = summary_proj, mapping = aes(x = year, y = median, color = 'Projection')) + 
    geom_line(data = summary_proj, mapping = aes(x = year, y = median, color = 'Projection')) + 
    geom_line(data = summary_proj, mapping = aes(x = year, y = median, color = 'Projection')) + 
    geom_line(data = summary_proj, mapping = aes(x = year, y = median, color = 'Projection')) + 
    geom_line(data = subset(all_hist_df, model == model_x), mapping = aes(x = year, y = temperature, color = 'Model'), size = 2) + 
    scale_color_manual(values = c('HadCrut4' = 'black', 'True' = 'red', 'Projection' = 'blue', 'Model' = 'purple')) + 
    geom_line(aes(y = observed_anomaly, color = 'HadCrut4')) + 
    geom_ribbon(aes(ymin = q025, ymax = q975), fill = 'red', alpha = 0.2) + 
    geom_ribbon(aes(ymin = q100, ymax = q900), fill = 'red', alpha = 0.3) + 
    geom_ribbon(aes(ymin = q025, ymax = q975), data = summary_proj,  fill = 'blue', alpha = 0.2) + 
    geom_ribbon(aes(ymin = q100, ymax = q900), data = summary_proj, fill = 'blue', alpha = 0.3) + 
    xlab('year') + ylab('Anomaly') + theme(text=element_text(size=17))
  pdf(file = paste0('CMIP5_model_forecasts/', model_x, '.pdf'), width = 8.5, height = 5)
  
  gridExtra::grid.arrange(p)
  dev.off()
}

##############
plot_hist_df <- all_hist_df %>% group_by(model) %>% mutate(new_temperature = temperature - mean(temperature[1:20]))

gg_color_hue <- function(n) {
  hues = seq(15, 375, length=n+1)
  hcl(h=hues, l=65, c=100)[1:n]
}

# Dynamically generate default color values, but have Paid="black".
adj_names = sort(setdiff(unique(subset(plot_hist_df, model %in% all_models[1:10])$model), "HadCrut4"))
values = gg_color_hue(length(adj_names))
names(values) = adj_names
values = c(values, c(HadCrut4="black"))

plot_cmip5 <-  ggplot(data = subset(plot_hist_df, model %in% all_models[1:10]), 
                      aes(x = year, y = new_temperature, color = model)) + geom_line() + ylab("Temperature") + 
  geom_line(data = real_data, aes(x=V1, y=V6-mean(real_data$V6[12:31]), color="HadCrut4")) + ggtitle('Past Data') + 
  scale_colour_manual(values=values) + 
  theme(text=element_text(size=17))
  
pacf_bias <- plot_hist_df %>% 
  group_by(model) %>% 
  summarise(lag1 = pacf(temperature - real_data)$acf[1,1,1], 
            lag2 = pacf(temperature - real_data)$acf[2,1,1],
            lag3 = pacf(temperature - real_data)$acf[3,1,1])


make_bugs_data <- function(all_hist_df)
{
  year <- 1861:2005
  model_forecast <- all_hist_df$temperature
  hadCrut_observation <- all_hist_df$real_data[1:length(year)]
  sd_hadCrut <- real_data$V3[real_data$V1 %in% year]
  n_models <- length(unique(all_hist_df$model))
  model_forecast_anomaly <- model_forecast
  observed_anomaly <- hadCrut_observation
  n_years <- length(year)
  
  return (list(n_years = n_years, 
               sd_delta = sd_hadCrut,
               model_forecast = model_forecast_anomaly,
               observed_anomaly = observed_anomaly, 
               n_models = n_models))
}

fit_ar1_model <- function(all_hist_df,
                          n.iterations = 11000, n.adapt = 1000, n.chains = 3, thin = 1, 
                          model_name = 'bayes_ar1_new.bug',
                          var.list = c("true_anomaly", "rho", "sd_w"))
{
  library(rjags)
  library(coda)
  library(reshape2)
  print("Calling make.bugs.data")
  jags.input <- make_bugs_data(all_hist_df)
  # browser()
  print("Successfully called make.bugs.data")
  jags.out <- jags.model(model_name, 
                         data=jags.input,
                         n.chains=n.chains, n.adapt=n.adapt)
  jags.output <- coda.samples(jags.out,
                              var.list,
                              n.iter=n.iterations,
                              thin=thin)
  return (list(input = jags.input, output = jags.output))
  
}

fit_results_test <- fit_ar1_model(all_hist_df)

true_anomaly_names <- paste('true_anomaly[', 1:145, ']', sep = '')
posterior_anomaly <- as.matrix(fit_results_test[[2]][[1]][, true_anomaly_names])
summary_anomaly <- data.frame(year = 1861:2005)
summary_anomaly$median <- apply(posterior_anomaly, 2, median)
summary_anomaly$q025 <- apply(posterior_anomaly, 2, quantile, probs = 0.025)
summary_anomaly$q100 <- apply(posterior_anomaly, 2, quantile, probs = 0.1)
summary_anomaly$q900 <- apply(posterior_anomaly, 2, quantile, probs = 0.9)
summary_anomaly$q975 <- apply(posterior_anomaly, 2, quantile, probs = 0.975)
summary_anomaly$observed_anomaly <- fit_results_test$input$observed_anomaly

ggplot(data = summary_anomaly, aes(x = year))  +
  geom_line(aes(y = median, color = 'True'), size = 2) + scale_color_manual(values = c('HadCrut4' = 'black', 'True' = 'red')) + 
  geom_line(aes(y = observed_anomaly, color = 'HadCrut4')) + 
  geom_ribbon(aes(ymin = q025, ymax = q975), fill = 'red', alpha = 0.2) + 
  geom_ribbon(aes(ymin = q100, ymax = q900), fill = 'red', alpha = 0.3) + xlab('year') + ylab('Anomaly')




