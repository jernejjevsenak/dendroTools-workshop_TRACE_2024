library(dplR)
library(dendroTools)
library(dplyr)
library(ggpubr)
library(reshape2)

chronology <- read.rwl("Sez.rwl")
detrended <- chron(detrend(chronology, method = "Spline", nyrs = 33), prewhiten = TRUE)
residual <- detrended[,2, drop = FALSE]

precipitation <- read.csv("Psum_SEZ.txt", sep = " ", header = TRUE)
precipitation <- select(precipitation, date, p_sum)

daily_precipitation <- data_transform(precipitation, format = "daily")
monthly_precipitation <- data_transform(precipitation, format = "monthly")

temperature <- read.csv("Tavg_SEZ.txt", sep = " ", header = TRUE)
temperature <- select(temperature, date, t_avg)

daily_temperature <- data_transform(temperature, format = "daily")
monthly_temperature <- data_transform(temperature, format = "monthly")

glimpse_daily_data(daily_temperature)

dsp_daily <- daily_response_seascorr(response = residual, 
                            env_data_primary = daily_precipitation,
                            env_data_control = daily_temperature,
                            row_names_subset = TRUE, remove_insignificant = FALSE,
                            alpha = 0.05)

dsp_monthly <- monthly_response_seascorr(response = residual, 
                           env_data_primary = monthly_precipitation,
                           env_data_control = monthly_temperature,
                           row_names_subset = TRUE, remove_insignificant = FALSE,
                           alpha = 0.05)

p1 <- plot(dsp_daily, type = 2)
p2 <- plot(dsp_monthly, type = 2)

ggarrange(p1, p2, ncol = 1)
ggsave("compare_daily_monthly_partial.png", width = 20, height = 15)


# let's use the summary() to compare the optimal windows
summary(dsp_daily)
summary(dsp_monthly)


dsp_monthly$optimized_return



################################################################################
# Compare partial vs. classical correlation coefficients

ds_daily <- daily_response(response = residual, 
                                     env_data = daily_precipitation,
                                     row_names_subset = TRUE, remove_insignificant = FALSE,
                                     alpha = 0.05, previous_year = TRUE)

ds_monthly <- monthly_response(response = residual, 
                                         env_data = monthly_precipitation,
                                         row_names_subset = TRUE, remove_insignificant = FALSE,
                                         alpha = 0.05, previous_year = TRUE)

p3 <- plot(ds_daily, type = 2)
p4 <- plot(ds_monthly, type = 2)

ggarrange(p1, p3, ncol = 1)
ggsave("compare_partial_classical_daily.png", width = 20, height = 15)

ggarrange(p2, p4, ncol = 1)
ggsave("compare_partial_classical_monthly.png", width = 20, height = 15)







#############################################################
# Subset the analyses using day_interval and month_interval #
#############################################################

dss_daily <- daily_response(response = residual, 
                           env_data = daily_precipitation,
                           day_interval = c(-150, 300), previous_year = TRUE,
                           row_names_subset = TRUE, remove_insignificant = FALSE,
                           alpha = 0.05)

dss_monthly <- monthly_response(response = residual, 
                               env_data = monthly_precipitation,
                               month_interval = c(-5, 10), previous_year = TRUE,
                               row_names_subset = TRUE, remove_insignificant = FALSE,
                               alpha = 0.05)

?daily_response_seascorr
dss_monthly$optimized_return

p5 <- plot(dss_daily, type = 2)
p6 <- plot(dss_monthly, type = 2)

ggarrange(p3, p5, ncol = 1)
ggsave("compare_subset_nosubset_classical_daily.png", width = 20, height = 15)

ggarrange(p4, p6, ncol = 1)
ggsave("compare_subset_nosubset_monthly.png", width = 20, height = 15)

