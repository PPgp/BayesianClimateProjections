library(ncdf4)
library(ggplot2)
library(pracma)

# load data
setwd("C:/Users/xchen/Desktop/Study/Research/Adrian Raftery/LiuRaftery2020Code/CMIP5data")
CMIP5_model = list()
CMIP5_data = array(0, c(144, 72, 240, 39))
for (i in 1:10) {
  CMIP5_model[[i]] = nc_open(paste(c("cmip5_tas_Amon_mod_rcp85_1_mean_00", i-1, ".nc"), collapse = ""))
  CMIP5_data[,,,i] = ncvar_get(CMIP5_model[[i]], varid = "tas")
}
for (i  in 11:38) {
  CMIP5_model[[i]] = nc_open(paste(c("cmip5_tas_Amon_mod_rcp85_1_mean_0", i-1, ".nc"), collapse = ""))
  CMIP5_data[,,,i] = ncvar_get(CMIP5_model[[i]], varid = "tas")
}
# deal with model 39
CMIP5_model[[39]] = nc_open("cmip5_tas_Amon_mod_rcp85_1_mean_038.nc")
model39_data_monthly = ncvar_get(CMIP5_model[[39]], varid = "tas")
model39_data_yearly = array(0, c(144,72,240))
for (i in 1:240) {
  tmp = model39_data_monthly[,, (12*(i-1)+1): (12*i)]
  model39_data_yearly[,,i] = apply(tmp, c(1,2), mean)
}
CMIP5_data[,,,39] = model39_data_yearly
for (i in 1:39) {
  print(ncatt_get(CMIP5_model[[i]], varid = 0, attname = "model_id")$value)
}

CMIP5_data_adjusted = array(0, c(72, 144, 240, 39))
for (t in 1:240) {
  for (k in 1:39) {
    CMIP5_data_adjusted[,,t,k] = t(CMIP5_data[,,t,k])[72:1,]
  }
}

# local temp
setwd("C:/Users/xchen/Desktop/Study/Research/Adrian Raftery/LiuRaftery2020Code")
data = read.csv("era5data2020_250.csv", header = T)
data = data[,-1]
data.array = array(0, c(72,144,41))
for (i in 1:72) {
  for (j in 1:144) {
    for (l in 1:41) {
      data.array[i,j,l] = data[(l-1)*72+i,j]
      index = (i-1)*144+j
    }
  }
}

temp.2019 = matrix(0, 72, 145)
temp.2019[,1:144] = data.array[,,41] - 273.15
temp.2019[,145] = temp.2019[,1]
write.csv(temp.2019, file = "temp2019.csv")


# error: detrend
start.year = 1979
end.year = 2015
error = array(0, c(72, 145, end.year-start.year+1))
for (i in 1:72) {
  for (j in 1:144) {
    error[i,j,] = detrend(data.array[i,j,1:37] - mean(data.array[i,j,2:21]))
  }
}
acf_all = matrix(0, 72, 144)
for (i in 1:72) {
  for (j in 1:144) {
    acf_all[i,j] = acf(error[i,j,], plot = F)$acf[3,1,1]
  }
}

sd_all = matrix(0, 72, 145)
for (i in 1:72) {
  for (j in 1:144) {
    sd_all[i,j] = sd(error[i,j,])
  }
}
sd_all[,145] = sd_all[,1]

par(pin = c(5,3))
plot(1979:2015, data.array[20,20,1:37]-273.15,type = "l", xlab = "Year", ylab = expression(paste("Temperature ", "[ "^o, "C", "]")), yaxt='n')
axis(2, c(10,11,12,13))
lines(1979:2015,data.array[20,20,1:37] - error[20,20,]-273.15, type = "l", col=2)
points(2015, data.array[20,20,37]-273.15, pch=19)
text(2015, data.array[20,20,37]+0.2-273.15, expression(T[A]))
points(2015, data.array[20,20,37] - error[20,20,37]-273.15, pch=19, col=2)
text(2015, data.array[20,20,37]- error[20,20,37]-0.2-273.15, expression(T[B]), col=2)

hist(c(acf_all), main = "acf(2)", xlab = NA)


write.csv(sd_all, "sd_new_new.csv")
# error:loess
tmp_x = 1:37
tmp_error = array(0, c(72, 145, end.year-start.year+1))
for (i in 1:72) {
  for (j in 1:144) {
    tmp_error[i,j,] = loess(data.array[i,j,1:37]~tmp_x)$residuals
  }
}
plot(1979:2015, data.array[20,20,1:37]-273.15,type = "l", xlab = "Year", ylab = expression(paste("Temperature ", "[ "^o, "C", "]")), yaxt='n')
axis(2, c(10,11,12,13))
lines(1979:2015,data.array[20,20,1:37] - tmp_error[20,20,]-273.15, type = "l", col=2)

tmp_acf_all = matrix(0, 72, 144)
for (i in 1:72) {
  for (j in 1:144) {
    acf_all[i,j] = pacf(tmp_error[i,j,], plot = F)$acf[2,1,1]
  }
}
tmp_sd_all = matrix(0, 72, 145)
for (i in 1:72) {
  for (j in 1:144) {
    tmp_sd_all[i,j] = sd(tmp_error[i,j,])
  }
}
tmp_sd_all[,145] = tmp_sd_all[,1]
write.csv(tmp_sd_all, "tmp_sd_new_new.csv")


# projections

# we change 39 cmip5 models to 35 models

CMIP5_data_adjusted_old = CMIP5_data_adjusted
CMIP5_data_adjusted = CMIP5_data_adjusted_old[,,,-c(20,21,23,24)]

######################################################################################
# verify natural variability 
# not run for main codes
start.year = 2064
end.year = 2100
error = array(0, c(72, 145, end.year-start.year+1))
for (i in 1:72) {
  for (j in 1:144) {
    error[i,j,] = detrend(CMIP5_data_adjusted[i,j,204:240,20])
  }
}

sd_all = matrix(0, 72, 145)
for (i in 1:72) {
  for (j in 1:144) {
    sd_all[i,j] = sd(error[i,j,])
  }
}
sd_all[,145] = sd_all[,1]
write.csv(sd_all, "sd_new_new_verify.csv")

###################################################################################################



local_alltime_mean = apply(CMIP5_data_adjusted, c(1,2,3),mean)
times = 1861:2100
local_minus = matrix(0, 72,145)
local_minus[,1:144] = apply(local_alltime_mean[,,which(times==2080):which(times==2099)],c(1,2), mean)-apply(local_alltime_mean[,,which(times==1980):which(times==1999)],c(1,2), mean)
local_minus[,145] = local_minus[,1]
write.csv(local_minus, "local_minus_new_new_35.csv")
write.csv(local_minus, "local_minus_new_new.csv")

#######################################################################################
# change to 35 models

local_minus_model = array(0, c(72,145,35))
for (k in 1:35) {
  local_minus_model[,1:144,k] = apply(CMIP5_data_adjusted[,,which(times==2080):which(times==2099),k],c(1,2), mean)-apply(CMIP5_data_adjusted[,,which(times==1980):which(times==1999),k],c(1,2), mean)
}
local_minus_model[,145,] = local_minus_model[,1,]

lat = seq(88.75, -88.75,  -2.5)
lon = seq(1.25,358.75,2.5)
weight = cos(lat*pi/180)/sum(cos(lat*pi/180))
global_temp = rep(0, 240)
for (t in 1:240) {
  global_temp[t] = sum(local_alltime_mean[,,t]*weight)/144
}
global_minus = mean(global_temp[which(times==2080):which(times==2099)]) - mean(global_temp[which(times==1980):which(times==1999)])
pattern_new = local_minus/global_minus
write.csv(pattern_new, "pattern_new_new_35.csv")

global_temp_model = matrix(0, 240, 35)
global_minus_model = rep(0,35)
pattern_new_model = array(0, c(72,145,35))
for (k in 1:35) {
  for (t in 1:240) {
    global_temp_model[t,k] = sum(CMIP5_data_adjusted[,,t,k]*weight)/144
  }
  global_minus_model[k] = mean(global_temp_model[which(times==2080):which(times==2099),k]) - mean(global_temp_model[which(times==1980):which(times==1999),k])
  pattern_new_model[,,k] = local_minus_model[,,k]/global_minus_model[k]
}

#############################################################################################

local_minus_model = array(0, c(72,145,39))
for (k in 1:39) {
  local_minus_model[,1:144,k] = apply(CMIP5_data_adjusted[,,which(times==2080):which(times==2099),k],c(1,2), mean)-apply(CMIP5_data_adjusted[,,which(times==1980):which(times==1999),k],c(1,2), mean)
}
local_minus_model[,145,] = local_minus_model[,1,]

lat = seq(88.75, -88.75,  -2.5)
lon = seq(1.25,358.75,2.5)
weight = cos(lat*pi/180)/sum(cos(lat*pi/180))
global_temp = rep(0, 240)
for (t in 1:240) {
  global_temp[t] = sum(local_alltime_mean[,,t]*weight)/144
}
global_minus = mean(global_temp[which(times==2080):which(times==2099)]) - mean(global_temp[which(times==1980):which(times==1999)])
pattern_new = local_minus/global_minus
write.csv(pattern_new, "pattern_new_new.csv")

global_temp_model = matrix(0, 240, 39)
global_minus_model = rep(0,39)
pattern_new_model = array(0, c(72,145,39))
for (k in 1:39) {
  for (t in 1:240) {
    global_temp_model[t,k] = sum(CMIP5_data_adjusted[,,t,k]*weight)/144
  }
  global_minus_model[k] = mean(global_temp_model[which(times==2080):which(times==2099),k]) - mean(global_temp_model[which(times==1980):which(times==1999),k])
  pattern_new_model[,,k] = local_minus_model[,,k]/global_minus_model[k]
}

#load("projections_sampled.Rdata")
load("projections_sampled_model.Rdata")

#########################################################################################
# change 39 models to 35 models
projections_sampled_model_old = projections_sampled_model
projections_sampled_model = projections_sampled_model_old[,-c(3001:4000, (20*1000+1):(21*1000), (22*1000+1):(23*1000), (33*1000+1):(34*1000))]


load("real_data_adjusted.Rdata")
load("real_data.Rda")
### load real_data
#projections_sampled_adjusted = projections_sampled_model + mean(real_data$V6[132:156]) - mean(real_data$V6[12:31]) 
projections_sampled_adjusted = projections_sampled_model + mean(real_data$V6[132:156])
index_1000_model = c()
for (i in 1:35) {
  index_1000_model = c(index_1000_model, ((i-1)*1000+1):((i-1)*1000+200))
}
projections_sampled_adjusted = projections_sampled_adjusted[, index_1000_model]

# a_2015 = t(rep(real_data_adjusted[166,'V6'],1000))
# projections_sampled_adjusted = rbind(a_2015, projections_sampled_adjusted)
# global_mean.diff = diff(projections_sampled_adjusted)

local_temp_proj_quantile = array(0, c(85, 8, 72, 144))
local_temp_proj_quantile_model = array(0, c(85, 8, 72, 144))
local_temp_proj_quantile_model_nosd = array(0, c(85, 8, 72, 144))
local_temp_2081_2100_proj_quantile = array(0, c(7, 72,145))
local_temp_2081_2100_proj_quantile_model = array(0, c(7, 72,145))

###############################################################################################################
# change 39 models to 35 models

for (i in 1:72) {
  cat( "i: ", i, "\n")
  for (j in 1:144) {
    local_temp_proj_combined = matrix(0, 85, 7000)
    local_temp_proj_combined_model = matrix(0, 85, 7000)
    local_temp_proj_combined_model_nosd = matrix(0, 85, 7000)
    local_temp_2081_2100_proj_combined = rep(0,7000)
    local_temp_2081_2100_proj_combined_model = rep(0, 7000)
    
    for (k in 1:35) {
      for (l in 1:200) {
        res = rnorm(85,0,sd_all[i,j])
        local_temp_proj_combined[,(k-1)*200+l] = data.array[i,j,37]- error[i,j,37] + res+ 
          pattern_new[i,j]*(projections_sampled_adjusted[,(k-1)*200+l] - real_data$V6[166])
        #- error[i,j,37]
        
        local_temp_proj_combined_model[,(k-1)*200+l] = data.array[i,j,37] - error[i,j,37]+ res+ 
          pattern_new_model[i,j,k]*(projections_sampled_adjusted[,(k-1)*200+l] - real_data$V6[166])
        #- error[i,j,37]
        
        local_temp_proj_combined_model_nosd[, (k-1)*200+l] = data.array[i,j,37] - error[i,j,37]+ 
          pattern_new_model[i,j,k]*(projections_sampled_adjusted[,(k-1)*200+l] - real_data$V6[166])
        
      }
    }
    local_temp_2081_2100_proj_combined = apply(local_temp_proj_combined[66:85,], 2, mean)
    local_temp_2081_2100_proj_combined_model = apply(local_temp_proj_combined_model[66:85,], 2, mean)
    
    local_temp_2081_2100_proj_quantile[,i,j] = quantile(local_temp_2081_2100_proj_combined, probs=c(0.025,0.05, 0.1, 0.5, 0.9, 0.95,0.975))
    local_temp_2081_2100_proj_quantile_model[,i,j] = quantile(local_temp_2081_2100_proj_combined_model, probs=c(0.025, 0.05, 0.1, 0.5, 0.9, 0.95, 0.975))
    
    local_temp_proj_quantile[, 1:7, i,j] = t(apply(local_temp_proj_combined, 1, quantile, probs=c(0.025, 0.05, 0.1, 0.5, 0.9, 0.95, 0.975)))
    local_temp_proj_quantile[, 8, i,j] = 2016:2100
    
    local_temp_proj_quantile_model[, 1:7, i,j] = t(apply(local_temp_proj_combined_model, 1, quantile, probs=c(0.025, 0.05, 0.1, 0.5, 0.9, 0.95, 0.975)))
    local_temp_proj_quantile_model[, 8, i,j] = 2016:2100
    
    local_temp_proj_quantile_model_nosd[, 1:7, i,j] = t(apply(local_temp_proj_combined_model_nosd, 1, quantile, probs=c(0.025, 0.05, 0.1, 0.5, 0.9, 0.95, 0.975)))
    local_temp_proj_quantile_model_nosd[, 8, i,j] = 2016:2100
    
  }
}

save(local_temp_proj_quantile,local_temp_proj_quantile_model,local_temp_2081_2100_proj_quantile,local_temp_2081_2100_proj_quantile_model , file = "proj_nonadjusted_35.Rdata")
load("proj_nonadjusted_35.Rdata")
save(local_temp_proj_quantile_model_nosd, file="proj_nonadjusted_nosd_35.Rdata")
load("proj_nonadjusted_nosd_35.Rdata")

####################################################################################################################################

for (i in 1:72) {
  cat( "i: ", i, "\n")
  for (j in 1:144) {
    local_temp_proj_combined = matrix(0, 85, 7800)
    local_temp_proj_combined_model = matrix(0, 85, 7800)
    local_temp_proj_combined_model_nosd = matrix(0, 85, 7800)
    local_temp_2081_2100_proj_combined = rep(0,7800)
    local_temp_2081_2100_proj_combined_model = rep(0, 7800)
    
    for (k in 1:39) {
      for (l in 1:200) {
        res = rnorm(85,0,sd_all[i,j])
        local_temp_proj_combined[,(k-1)*200+l] = data.array[i,j,37]- error[i,j,37] + res+ 
          pattern_new[i,j]*(projections_sampled_adjusted[,(k-1)*200+l] - real_data$V6[166])
        #- error[i,j,37]
       
        local_temp_proj_combined_model[,(k-1)*200+l] = data.array[i,j,37] - error[i,j,37]+ res+ 
          pattern_new_model[i,j,k]*(projections_sampled_adjusted[,(k-1)*200+l] - real_data$V6[166])
        #- error[i,j,37]
        
        local_temp_proj_combined_model_nosd[, (k-1)*200+l] = data.array[i,j,37] - error[i,j,37]+ 
          pattern_new_model[i,j,k]*(projections_sampled_adjusted[,(k-1)*200+l] - real_data$V6[166])
        
      }
    }
    local_temp_2081_2100_proj_combined = apply(local_temp_proj_combined[66:85,], 2, mean)
    local_temp_2081_2100_proj_combined_model = apply(local_temp_proj_combined_model[66:85,], 2, mean)
    
    local_temp_2081_2100_proj_quantile[,i,j] = quantile(local_temp_2081_2100_proj_combined, probs=c(0.025,0.05, 0.1, 0.5, 0.9, 0.95,0.975))
    local_temp_2081_2100_proj_quantile_model[,i,j] = quantile(local_temp_2081_2100_proj_combined_model, probs=c(0.025, 0.05, 0.1, 0.5, 0.9, 0.95, 0.975))
    
    local_temp_proj_quantile[, 1:7, i,j] = t(apply(local_temp_proj_combined, 1, quantile, probs=c(0.025, 0.05, 0.1, 0.5, 0.9, 0.95, 0.975)))
    local_temp_proj_quantile[, 8, i,j] = 2016:2100
    
    local_temp_proj_quantile_model[, 1:7, i,j] = t(apply(local_temp_proj_combined_model, 1, quantile, probs=c(0.025, 0.05, 0.1, 0.5, 0.9, 0.95, 0.975)))
    local_temp_proj_quantile_model[, 8, i,j] = 2016:2100
    
    local_temp_proj_quantile_model_nosd[, 1:7, i,j] = t(apply(local_temp_proj_combined_model_nosd, 1, quantile, probs=c(0.025, 0.05, 0.1, 0.5, 0.9, 0.95, 0.975)))
    local_temp_proj_quantile_model_nosd[, 8, i,j] = 2016:2100
    
  }
}

save(local_temp_proj_quantile,local_temp_proj_quantile_model,local_temp_2081_2100_proj_quantile,local_temp_2081_2100_proj_quantile_model , file = "proj_nonadjusted.Rdata")
load("proj_nonadjusted.Rdata")
save(local_temp_proj_quantile_model_nosd, file="proj_nonadjusted_nosd.Rdata")
load("proj_nonadjusted_nosd.Rdata")

################################################################################
################# Compute the "projected change WITHOUT natural variability" 
################# divided by "the amplitude of natural variability (1 sigma of epsilon)"
baseline = apply(data.array[,,1:20], c(1,2), mean)

ProjectedChangeWithoutNV_05 = matrix(0,72,145)
ProjectedChangeWithoutNV_95 = matrix(0,72,145)
ProjectedChangeWithoutNV_05[,1:144] = local_temp_proj_quantile_model_nosd[85,2,,] - baseline
ProjectedChangeWithoutNV_05[,145] = ProjectedChangeWithoutNV_05[,1]
SNR_05 = ProjectedChangeWithoutNV_05/sd_all
write.csv(SNR_05, "SNR_05_35.csv")
write.csv(ProjectedChangeWithoutNV_05, "ProjectedChangeWithoutNV_05_35.csv")
write.csv(SNR_05, "SNR_05.csv")
write.csv(ProjectedChangeWithoutNV_05, "ProjectedChangeWithoutNV_05.csv")

ProjectedChangeWithoutNV_95[,1:144] = local_temp_proj_quantile_model_nosd[85,6,,] - baseline
ProjectedChangeWithoutNV_95[,145] = ProjectedChangeWithoutNV_95[,1]
SNR_95 = ProjectedChangeWithoutNV_95/sd_all
write.csv(SNR_95, "SNR_95_35.csv")
write.csv(ProjectedChangeWithoutNV_95, "ProjectedChangeWithoutNV_95_35.csv")
write.csv(SNR_95, "SNR_95.csv")
write.csv(ProjectedChangeWithoutNV_95, "ProjectedChangeWithoutNV_95.csv")

################################################################################

which(which(local_temp_proj_quantile_model[85,4,,] - 273.15 < 0, arr.ind = T)[,1] == 15)
which(local_temp_proj_quantile_model[85,4,,] - 273.15 < 0, arr.ind = T)[800,]
# for (i in 1:72) {
#   for (j in 1:144) {
#     if(j %% 10 ==0){print(j)} 
#     for (l in 1:39000) {
#       local_temp_proj_combined[,l,i,j] = mean(local_alltime_mean[i,j,which(times==1980):which(times==1999)])+ rnorm(85,0,sd_all[i,j])+ 
#         pattern_new[i,j]*(projections_sampled_adjusted[,l]+ 273.15 - mean(global_temp[which(times==1980):which(times==1999)]))
#     }
#   }
# }
# local_temp_proj_quantile = array(0, c(86, 6, 72,144))
# for (i in 1:72) {
#   for (j in 1:144) {
#     local_temp_proj_quantile[, 1:5, i,j] = t(apply(local_temp_proj_combined[,,i,j], 1, quantile, probs=c(0.05, 0.1, 0.5, 0.9, 0.95)))
#     local_temp_proj_quantile[, 6, i,j] = 2015:2100
#   }
# }


temp2100_new_new = matrix(0, 72,145)
temp2100_new_new[,1:144] = local_temp_proj_quantile[85,4,,] - 273.15
temp2100_new_new[,145] = temp2100_new_new[,1]
write.csv(temp2100_new_new, "temp2100_new_new_35.csv")
write.csv(temp2100_new_new, "temp2100_new_new.csv")

temp2100_new_new_model = matrix(0, 72,145)
temp2100_new_new_model[,1:144] = local_temp_proj_quantile_model[85,4,,] - 273.15
temp2100_new_new_model[,145] = temp2100_new_new_model[,1]
write.csv(temp2100_new_new_model, "temp2100_new_new_model2_35.csv")
write.csv(temp2100_new_new_model, "temp2100_new_new_model2.csv")
xx = read.csv("temp2100_new_new_model.csv")
temp2100_new_new_model - xx[,-1]


xxx = read.csv("temp2100_new_new_model.csv")[,-1]
# anomaly
temp2100_anomaly_new_new = matrix(0, 72, 145)
for (i in 1:72) {
  for (j in 1:144) {
    temp2100_anomaly_new_new[i,j] = temp2100_new_new[i,j] + 273.15 - mean(data.array[i,j,1:20])
  }
}
temp2100_anomaly_new_new[,145] = temp2100_anomaly_new_new[,1]
write.csv(temp2100_anomaly_new_new, "temp2100_anomaly_new_new_35.csv")
write.csv(temp2100_anomaly_new_new, "temp2100_anomaly_new_new.csv")


temp2100_anomaly_new_new_model = matrix(0, 72, 145)
for (i in 1:72) {
  for (j in 1:144) {
    temp2100_anomaly_new_new_model[i,j] = temp2100_new_new_model[i,j] + 273.15 - mean(data.array[i,j,1:20])
  }
}
temp2100_anomaly_new_new_model[,145] = temp2100_anomaly_new_new_model[,1]
write.csv(temp2100_anomaly_new_new_model, "temp2100_anomaly_new_new_model2_35.csv")
write.csv(temp2100_anomaly_new_new_model, "temp2100_anomaly_new_new_model2.csv")
xx = read.csv("temp2100_anomaly_new_new_model.csv")
temp2100_anomaly_new_new_model - xx[,-1]


which(temp2100_anomaly_new_new == min(temp2100_anomaly_new_new[,1:144]), arr.ind = T)
which(temp2100_anomaly_new_new < -2 , arr.ind = T)

temp2100_95_anomaly_new_new = matrix(0, 72, 145)
temp2100_05_anomaly_new_new = matrix(0, 72, 145)
for (i in 1:72) {
  for (j in 1:144) {
    temp2100_95_anomaly_new_new[i,j] = local_temp_proj_quantile[85,6,i,j] - mean(data.array[i,j,1:20])
    temp2100_05_anomaly_new_new[i,j] = local_temp_proj_quantile[85,2,i,j] - mean(data.array[i,j,1:20])
  }
}
temp2100_05_anomaly_new_new[,145] = temp2100_05_anomaly_new_new[,1]
temp2100_95_anomaly_new_new[,145] = temp2100_95_anomaly_new_new[,1]
write.csv(temp2100_95_anomaly_new_new, "temp2100_95_anomaly_new_new_35.csv")
write.csv(temp2100_05_anomaly_new_new, "temp2100_05_anomaly_new_new_35.csv")
write.csv(temp2100_95_anomaly_new_new, "temp2100_95_anomaly_new_new.csv")
write.csv(temp2100_05_anomaly_new_new, "temp2100_05_anomaly_new_new.csv")

temp2100_95_anomaly_new_new_model = matrix(0, 72, 145)
temp2100_05_anomaly_new_new_model = matrix(0, 72, 145)
for (i in 1:72) {
  for (j in 1:144) {
    temp2100_95_anomaly_new_new_model[i,j] = local_temp_proj_quantile_model[85,6,i,j] - mean(data.array[i,j,1:20])
    temp2100_05_anomaly_new_new_model[i,j] = local_temp_proj_quantile_model[85,2,i,j] - mean(data.array[i,j,1:20])
  }
}
temp2100_05_anomaly_new_new_model[,145] = temp2100_05_anomaly_new_new_model[,1]
temp2100_95_anomaly_new_new_model[,145] = temp2100_95_anomaly_new_new_model[,1]
CI_model = temp2100_95_anomaly_new_new_model - temp2100_05_anomaly_new_new_model
write.csv(temp2100_95_anomaly_new_new_model, "temp2100_95_anomaly_new_new_model2_35.csv")
write.csv(temp2100_05_anomaly_new_new_model, "temp2100_05_anomaly_new_new_model2_35.csv")
write.csv(CI_model, "CI_model_35.csv")
write.csv(temp2100_95_anomaly_new_new_model, "temp2100_95_anomaly_new_new_model2.csv")
write.csv(temp2100_05_anomaly_new_new_model, "temp2100_05_anomaly_new_new_model2.csv")
write.csv(CI_model, "CI_model.csv")

xx = read.csv("temp2100_95_anomaly_new_new_model.csv")
xx[,-1] - temp2100_95_anomaly_new_new_model

xx = read.csv("CI_model.csv")
CI_model - xx[,-1]

# 2050 anomaly
temp2050_95_anomaly_new_new_model = matrix(0, 72, 145)
temp2050_05_anomaly_new_new_model = matrix(0, 72, 145)
for (i in 1:72) {
  for (j in 1:144) {
    temp2050_95_anomaly_new_new_model[i,j] = local_temp_proj_quantile_model[35,6,i,j] - mean(data.array[i,j,1:20])
    temp2050_05_anomaly_new_new_model[i,j] = local_temp_proj_quantile_model[35,2,i,j] - mean(data.array[i,j,1:20])
  }
}
temp2050_05_anomaly_new_new_model[,145] = temp2050_05_anomaly_new_new_model[,1]
temp2050_95_anomaly_new_new_model[,145] = temp2050_95_anomaly_new_new_model[,1]
write.csv(temp2050_05_anomaly_new_new_model/sd_all, "temp2050_standard_35.csv")
write.csv(temp2100_05_anomaly_new_new_model/sd_all, "temp2100_standard_35.csv")
write.csv(temp2050_95_anomaly_new_new_model, "temp2050_95_anomaly_new_new_model_35.csv")
write.csv(temp2050_05_anomaly_new_new_model, "temp2050_05_anomaly_new_new_model_35.csv")

write.csv(temp2050_05_anomaly_new_new_model/sd_all, "temp2050_standard.csv")
write.csv(temp2100_05_anomaly_new_new_model/sd_all, "temp2100_standard.csv")
write.csv(temp2050_95_anomaly_new_new_model, "temp2050_95_anomaly_new_new_model.csv")
write.csv(temp2050_05_anomaly_new_new_model, "temp2050_05_anomaly_new_new_model.csv")

xx = read.csv("temp2100_standard.csv")
temp2100_05_anomaly_new_new_model/sd_all - xx[,-1]

# plot single grid box
plot_grid = function(i,j){
  data.plot = matrix(0, 86, 6)
  data.plot[2:86,] = local_temp_proj_quantile[,,i,j]
  data.plot[2:86,1:5] = data.plot[2:86,1:5] - mean(data.array[i,j,1:20])
  colnames(data.plot) = c(paste0('Quant', c(0.05, 0.1, 0.5, 0.9, 0.95)), "year")
  data.plot = data.frame(data.plot)
  data.hist = matrix(0,37,2)
  data.hist[,1] = 1979:2015
  data.hist[,2] = data.array[i,j,1:37] - mean(data.array[i,j,1:20])
  colnames(data.hist) = c("year","hist")
  data.hist = data.frame(data.hist)
  data.plot[1,6] = 2015
  data.plot[1,1:5] = rep(data.hist[37,2],5)
  
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
    geom_ribbon(aes(x=year, ymin=Quant0.05, ymax=Quant0.95), data=data.plot,
                alpha=0.2, fill=color.us) +
    geom_line(data=data.plot, aes(x=year, y=Quant0.5),
              alpha=1, size=1.3, color=color.us) +
    # scale_x_continuous(breaks=c(2010, seq(2025, year.end, 25)), expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0))
  
  plot.obj <- plot.obj + 
    geom_line(aes(x = year, y = hist), data = data.hist)
  plot.obj
}

plot_grid_model = function(i,j){
  data.plot_model = matrix(0, 85, 8)
  data.plot_model[1:85,] = local_temp_proj_quantile_model[,,i,j]
  data.plot_model[1:85,1:7] = data.plot_model[1:85,1:7] - mean(data.array[i,j,1:20])
  colnames(data.plot_model) = c(paste0('Quant', c(0.025,0.05, 0.1, 0.5, 0.9, 0.95,0.975)), "year")
  data.plot_model = data.frame(data.plot_model)
  data.hist = matrix(0,41,2)
  data.hist[,1] = 1979:2019
  data.hist[,2] = data.array[i,j,1:41] - mean(data.array[i,j,1:20])
  colnames(data.hist) = c("year","hist")
  data.hist = data.frame(data.hist)
  
  data.cmip = matrix(0,122,2)
  data.cmip[,1] = 1979:2100
  data.cmip[,2] = local_alltime_mean[i,j, 119:240] - mean(data.array[i,j,1:20])
  colnames(data.cmip) = c("year","hist")
  data.cmip = data.frame(data.cmip)
  #data.plot_model[1,8] = 2015
  #data.plot_model[1,1:7] = rep(data.hist[37,2],7)
  
  data.plot = matrix(0, 85, 8)
  data.plot[1:85,] = local_temp_proj_quantile[,,i,j]
  data.plot[1:85,1:7] = data.plot[1:85,1:7] - mean(data.array[i,j,1:20])
  colnames(data.plot) = c(paste0('Quant', c(0.025,0.05, 0.1, 0.5, 0.9, 0.95,0.975)), "year")
  data.plot = data.frame(data.plot)
  #data.plot[1,8] = 2015
  #data.plot[1,1:7] = rep(data.hist[37,2],7)
  
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
    geom_ribbon(aes(x=year, ymin=Quant0.05, ymax=Quant0.95), data=data.plot_model,
                alpha=0.2, fill="blue") +
    geom_line(data=data.plot_model, aes(x=year, y=Quant0.5),
              alpha=1, size=1.3, color="blue") +
    geom_ribbon(aes(x=year, ymin=Quant0.05, ymax=Quant0.95), data=data.plot,
                alpha=0.2, fill="red") +
    geom_line(data=data.plot, aes(x=year, y=Quant0.5),
              alpha=1, size=1.3, color="red") + 
  # scale_x_continuous(breaks=c(2010, seq(2025, year.end, 25)), expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0))
  
  plot.obj <- plot.obj + 
    geom_line(aes(x = year, y = hist), data = data.hist)
  # plot.obj <- plot.obj+
  #   geom_line(aes(x = year, y = hist), data = data.cmip, color = "green")
  plot.obj <- plot.obj + geom_point(aes(x = year, y = hist),data = data.hist[37:41,], size = 1.5)
  plot.obj
}

plot_grid_model_new = function(i,j, lower, upper, name){
  data.plot_model = matrix(0, 85, 8)
  data.plot_model[1:85,] = local_temp_proj_quantile_model[,,i,j]
  data.plot_model[1:85,1:7] = data.plot_model[1:85,1:7] - mean(data.array[i,j,1:20])
  colnames(data.plot_model) = c(paste0('Quant', c(0.025,0.05, 0.1, 0.5, 0.9, 0.95,0.975)), "year")
  data.plot_model = data.frame(data.plot_model)
  data.hist = matrix(0,41,2)
  data.hist[,1] = 1979:2019
  data.hist[,2] = data.array[i,j,1:41] - mean(data.array[i,j,1:20])
  colnames(data.hist) = c("year","hist")
  data.hist = data.frame(data.hist)
  
  data.cmip = matrix(0,122,2)
  data.cmip[,1] = 1979:2100
  data.cmip[,2] = local_alltime_mean[i,j, 119:240] - mean(data.array[i,j,1:20])
  colnames(data.cmip) = c("year","hist")
  data.cmip = data.frame(data.cmip)
  #data.plot_model[1,8] = 2015
  #data.plot_model[1,1:7] = rep(data.hist[37,2],7)
  
  data.plot = matrix(0, 85, 8)
  data.plot[1:85,] = local_temp_proj_quantile[,,i,j]
  data.plot[1:85,1:7] = data.plot[1:85,1:7] - mean(data.array[i,j,1:20])
  colnames(data.plot) = c(paste0('Quant', c(0.025,0.05, 0.1, 0.5, 0.9, 0.95,0.975)), "year")
  data.plot = data.frame(data.plot)
  #data.plot[1,8] = 2015
  #data.plot[1,1:7] = rep(data.hist[37,2],7)
  
  font <- 'Times'
  plot.obj <- ggplot() +
    labs(color="Scenario") +
    theme(panel.background=element_rect(fill="white"),
          panel.grid.minor=element_line(color="white"),
          panel.grid.major=element_line(color="grey90"),
          text=element_text(family=font, size=22)) + ylim(lower, upper) + labs(title = name)
  
  color.us <- "#e41a1c"
  plot.obj <- plot.obj +
    ylab('Anomaly') +
    geom_ribbon(aes(x=year, ymin=Quant0.05, ymax=Quant0.95), data=data.plot_model,
                alpha=0.2, fill="blue") +
    geom_line(data=data.plot_model, aes(x=year, y=Quant0.5),
              alpha=1, size=1.3, color="blue") +
    geom_ribbon(aes(x=year, ymin=Quant0.05, ymax=Quant0.95), data=data.plot,
                alpha=0.2, fill="red") +
    geom_line(data=data.plot, aes(x=year, y=Quant0.5),
              alpha=1, size=1.3, color="red") 
    # scale_x_continuous(breaks=c(2010, seq(2025, year.end, 25)), expand = c(0,0)) +
    #scale_y_continuous(expand = c(0,0))
  
  plot.obj <- plot.obj + 
    geom_line(aes(x = year, y = hist), data = data.hist)
  # plot.obj <- plot.obj+
  #   geom_line(aes(x = year, y = hist), data = data.cmip, color = "green")
  plot.obj <- plot.obj + geom_point(aes(x = year, y = hist),data = data.hist[37:41,], size = 1.5)
  plot.obj
}

plot_grid_model_new_new = function(i,j, lower, upper, name){
  data.plot_model = matrix(0, 85, 8)
  data.plot_model[1:85,] = local_temp_proj_quantile_model[,,i,j]
  data.plot_model[1:85,1:7] = data.plot_model[1:85,1:7] - mean(data.array[i,j,1:20])
  colnames(data.plot_model) = c(paste0('Quant', c(0.025,0.05, 0.1, 0.5, 0.9, 0.95,0.975)), "year")
  data.plot_model = data.frame(data.plot_model)
  
  blue1 = data.plot_model[1,6] - data.plot_model[1,2]
  blue2 = data.plot_model[85,6] - data.plot_model[85,2]
  
  data.plot_model_nosd = matrix(0, 85, 8)
  data.plot_model_nosd[1:85,] = local_temp_proj_quantile_model_nosd[,,i,j]
  data.plot_model_nosd[1:85,1:7] = data.plot_model_nosd[1:85,1:7] - mean(data.array[i,j,1:20])
  colnames(data.plot_model_nosd) = c(paste0('Quant', c(0.025,0.05, 0.1, 0.5, 0.9, 0.95,0.975)), "year")
  data.plot_model_nosd = data.frame(data.plot_model_nosd)
  
  green1 = data.plot_model_nosd[1,6] - data.plot_model_nosd[1,2]
  green2 = data.plot_model_nosd[85,6] - data.plot_model_nosd[85,2]
  
  data.hist = matrix(0,41,2)
  data.hist[,1] = 1979:2019
  data.hist[,2] = data.array[i,j,1:41] - mean(data.array[i,j,1:20])
  colnames(data.hist) = c("year","hist")
  data.hist = data.frame(data.hist)
  
  data.cmip = matrix(0,122,2)
  data.cmip[,1] = 1979:2100
  data.cmip[,2] = local_alltime_mean[i,j, 119:240] - mean(data.array[i,j,1:20])
  colnames(data.cmip) = c("year","hist")
  data.cmip = data.frame(data.cmip)
  #data.plot_model[1,8] = 2015
  #data.plot_model[1,1:7] = rep(data.hist[37,2],7)
  
  data.plot = matrix(0, 85, 8)
  data.plot[1:85,] = local_temp_proj_quantile[,,i,j]
  data.plot[1:85,1:7] = data.plot[1:85,1:7] - mean(data.array[i,j,1:20])
  colnames(data.plot) = c(paste0('Quant', c(0.025,0.05, 0.1, 0.5, 0.9, 0.95,0.975)), "year")
  data.plot = data.frame(data.plot)
  #data.plot[1,8] = 2015
  #data.plot[1,1:7] = rep(data.hist[37,2],7)
  
  font <- 'Times'
  plot.obj <- ggplot() +
    labs(color="Scenario") +
    theme(panel.background=element_rect(fill="white"),
          panel.grid.minor=element_line(color="white"),
          panel.grid.major=element_line(color="grey90"),
          text=element_text(family=font, size=22)) + ylim(lower, upper) + labs(title = name)
  
  color.us <- "#e41a1c"
  plot.obj <- plot.obj +
    ylab('Anomaly') +
    geom_ribbon(aes(x=year, ymin=Quant0.05, ymax=Quant0.95), data=data.plot_model,
                alpha=0.2, fill="blue") +
    geom_line(data=data.plot_model, aes(x=year, y=Quant0.5),
              alpha=1, size=1.3, color="blue") +
    geom_ribbon(aes(x=year, ymin=Quant0.05, ymax=Quant0.95), data=data.plot,
                alpha=0.2, fill="red") +
    geom_line(data=data.plot, aes(x=year, y=Quant0.5),
              alpha=1, size=1.3, color="red") +
    geom_line(data=data.plot_model_nosd, aes(x=year, y=Quant0.95),
              alpha=1, size=1.3, color="green") + 
    geom_line(data=data.plot_model_nosd, aes(x=year, y=Quant0.05),
              alpha=1, size=1.3, color="green") 
    
  # scale_x_continuous(breaks=c(2010, seq(2025, year.end, 25)), expand = c(0,0)) +
  #scale_y_continuous(expand = c(0,0))
  
  plot.obj <- plot.obj + 
    geom_line(aes(x = year, y = hist), data = data.hist)
  # plot.obj <- plot.obj+
  #   geom_line(aes(x = year, y = hist), data = data.cmip, color = "green")
  plot.obj <- plot.obj + geom_point(aes(x = year, y = hist),data = data.hist[37:41,], size = 1.5)
  plot.obj <- plot.obj + geom_text(aes(x=2100,y=data.plot_model[85,4]+0.1,label=round(data.plot_model[85,4],1)))
  plot.obj <- plot.obj + geom_text(aes(x=2016,y=data.plot_model[1,4]+0.3,label=round(data.plot_model[1,4],1)))
  plot.obj <- plot.obj + geom_text(aes(x=2100,y=data.plot_model[85,2]+0.1,label=round(data.plot_model[85,2],1)))
  plot.obj <- plot.obj + geom_text(aes(x=2100,y=data.plot_model[85,6]+0.1,label=round(data.plot_model[85,6],1)))
  plot.obj

  #return(list(green1 = green1, green2 = green2, blue1 = blue1, blue2 = blue2))
}




which(sd_all == max(sd_all), arr.ind = T)
which(sd_all == min(sd_all[,1:144]), arr.ind = T)
plot_grid(11,122)
plot_grid(34,53)

plot_grid_model(11,122)
plot_grid_model(34,53)

which(local_interval_length == max(local_interval_length), arr.ind = T)
which(local_interval_length == min(local_interval_length[,1:144]), arr.ind = T)
plot_grid(5,29)
plot_grid(56,53)
plot_grid_model(5,29)
plot_grid_model(56,53)

plot_grid_model_new(17,121,-3,8, "Seattle") #seattle
plot_grid_model_new(5, 6,-2.5,12, "Longyearbyen") #northernmost longyearbyen
plot_grid_model_new(37,104,-1.5,4, "Quito") #equator quito
plot_grid_model_new(58, 99, -1.5,4, "Puerto Williams") #southernmost Puerto Williams
plot_grid_model_new(17,1, -3,8, "Paris") #paris
plot_grid_model_new(20,108, -3,8, "Chicago") #chicago
plot_grid_model_new(14,16, -2.5, 8, "Moscow") #moscow
plot_grid_model_new(17,1 ) #paris
plot_grid_model_new(38,7, -1.5,5.5, "Kinshasha") #kinshasha
plot_grid_model_new(25,31, -3,8, " New Delhi") #new delhi
plot_grid_model_new(38,120, -1.5,5.5, "Manaus") #manaus
plot_grid_model_new(21, 47, -2.5,8, "Beijing") #beijing

seattle = plot_grid_model_new_new(17,121,-3,8, "Seattle") #seattle
longyearbyen = plot_grid_model_new_new(5, 6,-2.5,12, "Longyearbyen") #northernmost longyearbyen
quito = plot_grid_model_new_new(37,104,-1.5,5.5, "Quito") #equator quito
puerto = plot_grid_model_new_new(58, 99, -1.5,5.5, "Puerto Williams") #southernmost Puerto Williams
paris = plot_grid_model_new_new(17,1, -3,8, "Paris") #paris
chicago = plot_grid_model_new_new(20,108, -3,8, "Chicago") #chicago
moscow = plot_grid_model_new_new(14,16, -2.5, 8, "Moscow") #moscowplot_grid_model_new_new(17,1 ) #paris
kinshasa = plot_grid_model_new_new(38,7, -1.5,5.5, "Kinshasa") #kinshasha
newdelhi = plot_grid_model_new_new(25,31, -3,8, " New Delhi") #new delhi
manaus = plot_grid_model_new_new(38,120, -1.5,5.5, "Manaus") #manaus
beijing = plot_grid_model_new_new(21, 47, -2.5,8, "Beijing") #beijing
plot_grid_model_new_new(25,111,-3,8, "Houston")
plot_grid_model_new_new(50,61,-1.5,5.5, "Sydney")
plot_grid_model_new_new(34,2,-1.5,5.5, "Lagos")
plot_grid_model_new_new(24,13,-1.5,5.5, "Cairo")
plot_grid_model_new_new(31,49, -1.5,5.5, "Manila") 


city1 = c("Beijing", "Seattle", "Paris", "Kinshasa","Quito")
city2 = c("Moscow", "Chicago", "New Delhi", "Manaus", "Puerto Williams")
ratio11 = c(beijing$blue2/beijing$blue1, seattle$blue2/seattle$blue1, 
            paris$blue2/paris$blue1, kinshasa$blue2/kinshasa$blue1, quito$blue2/quito$blue1)
ratio12 = c(beijing$green2/beijing$green1, seattle$green2/seattle$green1, 
            paris$green2/paris$green1, kinshasa$green2/kinshasa$green1, quito$green2/quito$green1)
ratio21 = c(moscow$blue2/moscow$blue1, chicago$blue2/chicago$blue1,
            newdelhi$blue2/newdelhi$blue1, manaus$blue2/manaus$blue1, puerto$blue2/puerto$blue1)
ratio22 = c(moscow$green2/moscow$green1, chicago$green2/chicago$green1,
            newdelhi$green2/newdelhi$green1, manaus$green2/manaus$green1, puerto$green2/puerto$green1)


dframe1 = data.frame(City = city1, Ratio1 = ratio11, Ratio2 = ratio12)
dframe2 = data.frame(City = city2, Ratio1 = ratio21, Ratio2 = ratio22)


plot_grid_model(17,49) #seattle
plot_grid_model(5, 6) #northernmost longyearbyen
plot_grid_model(37,104) #equator quito
plot_grid_model(58, 99) #southernmost Puerto Williams



local_interval_length = matrix(0, 72, 145)
for (i in 1:72) {
  for (j in 1:144) {
    local_interval_length[i,j] = local_temp_proj_quantile[85,5,i,j] - local_temp_proj_quantile[85,1,i,j]
  }
}
local_interval_length[,145] = local_interval_length[,1]
write.csv(local_interval_length, "local_interval_length.csv")

local_interval_length_model = matrix(0, 72, 145)
for (i in 1:72) {
  for (j in 1:144) {
    local_interval_length_model[i,j] = local_temp_proj_quantile_model[85,5,i,j] - local_temp_proj_quantile_model[85,1,i,j]
  }
}
local_interval_length_model[,145] = local_interval_length_model[,1]
ratio_model = (local_interval_length_model - local_interval_length)/local_interval_length_model
ratio_model2 = (local_interval_length_model - local_interval_length)/local_interval_length

write.csv(local_interval_length_model, "local_interval_length_model.csv")



# 20 year's mean
# local_temp_2081_2100_proj_combined = array(0, c(1000,72,144))
# for (i in 1:72) {
#   for (j in 1:144) {
#     for (l in 1:1000) {
#       local_temp_2081_2100_proj_combined[l,i,j] = mean(local_temp_proj_combined[67:86, l, i, j])
#     }
#   }
# }


local_temp_2081_2100_proj_quantile_anomaly = array(0, c(7, 72,145))
for (i in 1:72) {
  for (j in 1:144) {
    local_temp_2081_2100_proj_quantile_anomaly[, i,j] = local_temp_2081_2100_proj_quantile[,i,j] - mean(data.array[i,j,1:20])
  }
}
local_temp_2081_2100_proj_quantile_anomaly[,,145] = local_temp_2081_2100_proj_quantile_anomaly[,,1]
write.csv(local_temp_2081_2100_proj_quantile_anomaly[4,,], "temp_2081_2100_anomaly_new_new_35.csv")
write.csv(local_temp_2081_2100_proj_quantile_anomaly[2,,], "temp_2081_2100_05_anomaly_new_new_35.csv")
write.csv(local_temp_2081_2100_proj_quantile_anomaly[6,,], "temp_2081_2100_95_anomaly_new_new_35.csv")

write.csv(local_temp_2081_2100_proj_quantile_anomaly[4,,], "temp_2081_2100_anomaly_new_new.csv")
write.csv(local_temp_2081_2100_proj_quantile_anomaly[2,,], "temp_2081_2100_05_anomaly_new_new.csv")
write.csv(local_temp_2081_2100_proj_quantile_anomaly[6,,], "temp_2081_2100_95_anomaly_new_new.csv")

local_temp_2081_2100_proj_quantile_anomaly_model = array(0, c(7, 72,145))
for (i in 1:72) {
  for (j in 1:144) {
    local_temp_2081_2100_proj_quantile_anomaly_model[, i,j] = local_temp_2081_2100_proj_quantile_model[,i,j] - mean(data.array[i,j,1:20])
  }
}
local_temp_2081_2100_proj_quantile_anomaly_model[,,145] = local_temp_2081_2100_proj_quantile_anomaly_model[,,1]
write.csv(local_temp_2081_2100_proj_quantile_anomaly_model[4,,], "temp_2081_2100_anomaly_new_new_model_35.csv")
write.csv(local_temp_2081_2100_proj_quantile_anomaly_model[2,,], "temp_2081_2100_05_anomaly_new_new_model_35.csv")
write.csv(local_temp_2081_2100_proj_quantile_anomaly_model[6,,], "temp_2081_2100_95_anomaly_new_new_model_35.csv")

write.csv(local_temp_2081_2100_proj_quantile_anomaly_model[4,,], "temp_2081_2100_anomaly_new_new_model.csv")
write.csv(local_temp_2081_2100_proj_quantile_anomaly_model[2,,], "temp_2081_2100_05_anomaly_new_new_model.csv")
write.csv(local_temp_2081_2100_proj_quantile_anomaly_model[6,,], "temp_2081_2100_95_anomaly_new_new_model.csv")

######################## Change to natural variability #############################################################
xxx = read.csv("temp_2081_2100_anomaly_new_new_model_35.csv", header = T, )[,-1]
NaturalVariabilityToChange =  sd_all/xxx
NaturalVariabilityToChange_abs =  sd_all/abs(xxx)
write.csv(NaturalVariabilityToChange, "NaturalVariabilityToChange.csv")
write.csv(NaturalVariabilityToChange_abs, "NaturalVariabilityToChange_abs.csv")

ChangetoNaturalVariability = xxx/sd_all
write.csv(ChangetoNaturalVariability, "ChangetoNaturalVariability_35.csv")
write.csv(ChangetoNaturalVariability, "ChangetoNaturalVariability.csv")

####################################################################################################################
# 20 year's CI length
anomaly05 = read.csv("temp_2081_2100_05_anomaly_new_new_model_35.csv")[,-1]
anomaly95 = read.csv("temp_2081_2100_95_anomaly_new_new_model_35.csv")[,-1]
CIlength_20 = anomaly95 - anomaly05
write.csv(CIlength_20, "CIlength_20_35.csv")
write.csv(CIlength_20, "CIlength_20.csv")

local_temp_2081_2100_proj_quantile_anomaly_model[4,,] - xxx[,-1]
local_temp_2081_2100_proj_quantile_anomaly_model[2,,] - anomaly05
local_temp_2081_2100_proj_quantile_anomaly_model[6,,] - anomaly95
