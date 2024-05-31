library(dplR)
library(dendroTools)
library(dplyr)
library(ggpubr)
library(reshape2)
library(ggplot2)

chronology <- read.rwl("Sez.rwl")
detrended <- detrend(chronology, method = "Spline", nyrs = 33)

summary(detrended)

site_chron <- chron(detrended, prewhiten = TRUE)
residual <- site_chron[,"res", drop = FALSE]
residual <- residual[!is.na(residual$res), ,drop = FALSE]

precipitation <- read.csv("Psum_SEZ.txt", sep = " ", header = TRUE)
precipitation <- select(precipitation, date, p_sum)

daily_precipitation <- data_transform(precipitation, format = "daily")

glimpse_daily_data(daily_precipitation)

monthly_precipitation <- data_transform(precipitation, format = "monthly")


ds_daily <- daily_response(response = residual, 
                           env_data = daily_precipitation,
                           row_names_subset = TRUE, 
                           remove_insignificant = FALSE,
                           alpha = 0.05, 
                           day_interval = c(-150, 300),
                           lower_limit = 21, upper_limit = 180)


plot(ds_daily)

melted <- melt(data.frame(ds_daily$calculations))
melted$season_length <- seq(21, 180)

melted$variable <- as.numeric(gsub("X", "", melted$variable))
melted <- dplyr::filter(melted, !is.na(value))


ggplot(melted, aes(x = variable, y = season_length, fill = value)) + geom_tile() +
  scale_x_continuous(expand=c(0,0),breaks = c(300, 466, 566, 666), 
                     labels = c("300*", 100, 200, 300)) +
  scale_y_continuous(expand=c(0,0)) + theme_bw() +
  geom_vline(linetype = "dashed", xintercept = 366, alpha = 0.5) +
  scale_fill_gradient2() +
  xlab("Season Start") +
  ylab("Season Length") +
  ggtitle("1950 - 1980") + labs(fill = "cor")

ggsave("heatmap example.png", bg = "white")

plot(ds_daily, type = 2)
  














ds_monthly <- monthly_response(response = residual, 
                               env_data = monthly_precipitation,
                               month_interval = c(-5, 10),
                               row_names_subset = TRUE, remove_insignificant = FALSE)

class(ds_monthly)
class(ds_daily)

p1 <- plot(ds_daily, type = 2)
p2 <- plot(ds_monthly, type = 2)

ggarrange(p1, p2, ncol = 1)
ggsave("compare_daily_monthly.png", width = 20, height = 15)


# let's use the summary() to compare the optimal windows
summary(ds_daily)
summary(ds_monthly)




################################################################################
# The importance of reference_window

ds_monthly_start <- monthly_response(response = residual, env_data = monthly_precipitation,
                               row_names_subset = TRUE, remove_insignificant = FALSE,
                               reference_window = "start",
                               alpha = 0.05)

ds_monthly_end <- monthly_response(response = residual, env_data = monthly_precipitation,
                                     row_names_subset = TRUE, remove_insignificant = FALSE,
                                     reference_window = "end",
                                     alpha = 0.05)


plot(ds_monthly_start, type = 2)
ggsave("ms_start.png", width = 12, height = 7)


plot(ds_monthly_end, type = 2)
ggsave("ms_end.png", width = 12, height = 7)

################################################################################

ds_daily_start <- daily_response(response = residual, env_data = daily_precipitation,
                           row_names_subset = TRUE, remove_insignificant = FALSE,
                           alpha = 0.05, reference_window = "start")

ds_daily_middle <- daily_response(response = residual, env_data = daily_precipitation,
                           row_names_subset = TRUE, remove_insignificant = FALSE,
                           alpha = 0.05, reference_window = "middle")

ds_daily_end <- daily_response(response = residual, env_data = daily_precipitation,
                           row_names_subset = TRUE, remove_insignificant = FALSE,
                           alpha = 0.05, reference_window = "end")

plot(ds_daily_start, type = 2) + theme(plot.title = element_blank(), legend.position = "none")
ggsave("ds_start.png", width = 17, height = 7)

plot(ds_daily_middle, type = 2) + theme(plot.title = element_blank(), legend.position = "none")
ggsave("ds_middle.png", width = 17, height = 7)

plot(ds_daily_end, type = 2)+ theme(plot.title = element_blank(), legend.position = "none")
ggsave("ds_end.png", width = 17, height = 7)




##########################################
# Filtering non-significant correlations #
##########################################

ds_monthly_ns <- monthly_response(response = residual, env_data = monthly_precipitation,
                                     row_names_subset = TRUE, remove_insignificant = FALSE,
                                     reference_window = "start",
                                     alpha = 0.05)

ds_monthly_05 <- monthly_response(response = residual, env_data = monthly_precipitation,
                                  row_names_subset = TRUE, remove_insignificant = TRUE,
                                  reference_window = "start",
                                  alpha = 0.05)

ds_monthly_001 <- monthly_response(response = residual, env_data = monthly_precipitation,
                                  row_names_subset = TRUE, remove_insignificant = TRUE,
                                  reference_window = "start",
                                  alpha = 0.001)

ds_monthly_0001 <- monthly_response(response = residual, env_data = monthly_precipitation,
                                   row_names_subset = TRUE, remove_insignificant = TRUE,
                                   reference_window = "start",
                                   alpha = 0.0001)

p1 <- plot(ds_monthly_ns, type = 2) + theme(plot.title = element_blank(), legend.position = "none") + 
  scale_fill_gradient2(limits = c(-0.20, 0.60))

p2 <- plot(ds_monthly_05, type = 2) + theme(plot.title = element_blank(), legend.position = "none") + 
  scale_fill_gradient2(limits = c(-0.20, 0.60))

p3 <- plot(ds_monthly_001, type = 2) + theme(plot.title = element_blank(), legend.position = "none") + 
  scale_fill_gradient2(limits = c(-0.20, 0.60))

p4 <- plot(ds_monthly_0001, type = 2) + theme(plot.title = element_blank(), legend.position = "none") + 
  scale_fill_gradient2(limits = c(-0.20, 0.60))

ggarrange(p1,p2,p3,p4)
ggsave("ms_alpha.png", width = 12, height = 12)


##############################


ds_daily_ns <- daily_response(response = residual, env_data = daily_precipitation,
                                  row_names_subset = TRUE, remove_insignificant = FALSE,
                                  reference_window = "start",
                                  alpha = 0.05, previous_year = TRUE)


plot(ds_daily_ns, type = 2) +
    #scale_x_continuous(breaks = seq(100, 700, by = 100),
    #                label = c("100*", "200*", "300*", "34", "134", "234", "334")) 

   scale_x_continuous(breaks = c(100, 200, 300, 466, 566, 666),
                   label = c("100*", "200*", "300*", "100", "200", "300")) 

ds_daily_05 <- daily_response(response = residual, env_data = daily_precipitation,
                                  row_names_subset = TRUE, remove_insignificant = TRUE,
                                  reference_window = "start",
                                  alpha = 0.05)

ds_daily_001 <- daily_response(response = residual, env_data = daily_precipitation,
                                   row_names_subset = TRUE, remove_insignificant = TRUE,
                                   reference_window = "start",
                                   alpha = 0.001)

ds_daily_0001 <- daily_response(response = residual, env_data = daily_precipitation,
                                    row_names_subset = TRUE, remove_insignificant = TRUE,
                                    reference_window = "start",
                                    alpha = 0.0001)

p1 <- plot(ds_daily_ns, type = 2) + theme(plot.title = element_blank(), legend.position = "none") + 
  scale_fill_gradient2(limits = c(-0.25, 0.70))

p2 <- plot(ds_daily_05, type = 2) + theme(plot.title = element_blank(), legend.position = "none") + 
  scale_fill_gradient2(limits = c(-0.25, 0.70))

p3 <- plot(ds_daily_001, type = 2) + theme(plot.title = element_blank(), legend.position = "none") + 
  scale_fill_gradient2(limits = c(-0.25, 0.70))

p4 <- plot(ds_daily_0001, type = 2) + theme(plot.title = element_blank(), legend.position = "none") + 
  scale_fill_gradient2(limits = c(-0.25, 0.70))

ggarrange(p1,p2,p3,p4)
ggsave("ds_alpha.png", width = 12, height = 12)





##################################################
# Steps to manually create heatmaps using ggplot #
##################################################

# 1 Extract table with correlations
data_cor <- data.frame(ds_daily_ns$calculations)

# 2 Add information from on season_length which is stored in row.names
data_cor$season_length <- row.names(data_cor)

# 3 Transform season_length to numeric variable
data_cor$season_length <- as.numeric(data_cor$season_length)

# 4 Transform data from "wide" to "long" format using the melt function from reshape2
# use season_length as id variable
data_cor <- melt(data_cor, id.vars = "season_length")
head(data_cor)

# 5 Transform variable (season_start) to numeric variable
data_cor$season_start <- as.numeric(data_cor$variable)

# 6 Create heatmap using the ggplot2 syntax
ggplot(data_cor, aes(x = season_start, y = season_length, fill = value)) + geom_tile() 

# 7 Some fine tuning won't hurt
ggplot(dplyr::filter(data_cor,!is.na(value)), # remove NA's to get transparent background 
       aes(x = season_start, y = season_length, fill = value)) + geom_tile() +
  scale_fill_gradient2(low = "red", high = "blue") + 
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  theme_bw() +
  xlab("Day of Year") +
  ylab("Season Length") +
  theme(legend.title = element_blank(),
        legend.position = "bottom")

ggsave("ds_daily_manual_ggplot.pdf", width = 10, height = 7)

plot(ds_daily_ns, type = 2) + theme(legend.position = "none") +xlab("") + 
  scale_fill_gradient2(low = "red", high = "blue") + 
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  theme_bw() +
  xlab("Day of Year") +
  ylab("Season Length") +
  theme(legend.title = element_blank(),
        legend.position = "bottom")

# Repeat for monthly data as well

# 1 Extract table with correlations
data_cor <- data.frame(ds_monthly_ns$calculations)

# 2 Add information from on season_length which is stored in row.names
data_cor$season_length <- row.names(data_cor)

# 3 Transform season_length to numeric variable
data_cor$season_length <- as.numeric(data_cor$season_length)

# 4 Transform data from "wide" to "long" format using the melt function from reshape2
# use season_length as id variable
data_cor <- melt(data_cor, id.vars = "season_length")
head(data_cor)

# 5 Transform variable (season_start) to numeric variable
data_cor$season_start <- as.numeric(data_cor$variable)

# 6 Create heatmap using the ggplot2 syntax
ggplot(data_cor, aes(x = season_start, y = season_length, fill = value)) + geom_tile() 

# 7 Some fine tuning won't hurt
ggplot(dplyr::filter(data_cor,!is.na(value)), # remove NA's to get transparent background 
       aes(x = season_start, y = season_length, fill = value)) + geom_tile() +
  scale_fill_gradient2(low = "red", high = "blue") +
  scale_x_continuous(expand = c(0,0), breaks = seq(1,12), 
                     labels = c("J", "F", "M", "A", "M","J","J","A","S","O","N","D")) +
  scale_y_continuous(expand = c(0,0)) +
  theme_bw() +
  xlab("Month") +
  ylab("Season Length") +
  theme(legend.title = element_blank(),
        legend.position = "bottom")

ggsave("ds_monthly_manual_ggplot.png", width = 10, height = 7)





