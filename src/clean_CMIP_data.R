library(ggplot2)
library(reshape2)
library(plotly)

all_files <- list.files(pattern = '*.dat', path = 'CMIP_data')
string_part1 <- sapply(all_files, function(x) {return (strsplit(x, '_')[[1]][1])})
string_part2 <- sapply(all_files, function(x) {return (strsplit(x, '_')[[1]][2])})
string_part3 <- sapply(all_files, function(x) {return (strsplit(x, '_')[[1]][3])})
string_part4 <- sapply(all_files, function(x) {return (strsplit(x, '_')[[1]][4])})
string_part5 <- sapply(all_files, function(x) {return (strsplit(x, '_')[[1]][5])})
string_part6 <- sapply(all_files, function(x) {return (strsplit(x, '_')[[1]][6])})
all_files_split <- data.frame(v1 = string_part1, v2 = string_part2, v3 = string_part3,
                              v4 = string_part4, v5 = string_part5, v6 = string_part6)

temp <- read.table('had4_krig_annual_v2_0_0.txt')
real_data <- temp

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

save(cleaned_all_data, file = 'cleaned_all_data.Rda')
save(real_data, file = 'real_data.Rda')




