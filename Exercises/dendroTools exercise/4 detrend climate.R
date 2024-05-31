library(dendroTools)
library(dplyr)
library(reshape2)
library(ggplot2)
library(lubridate)
library(dplyr)
library(broom)
library(dplR)

# A) Example with mean temperature

climate <- read.csv("Tavg_SEZ.txt", sep = " ", header = TRUE)
climate <- select(climate, date, t_avg)
tmp_model <- lm(t_avg ~ seq(1:nrow(climate)), data = climate)
climate$t_avg_predicted <- predict(tmp_model)

climate$residual <- climate$t_avg - climate$t_avg_predicted
climate$detrended <- climate$residual / sd(climate$residual)
climate <- select(climate, date, t_avg, detrended)

climate <- melt(climate, id.vars = "date")
climate$date <- ymd(climate$date)


p1 <- ggplot(climate, aes(x = date, y = value,  colour = variable)) + geom_line() +
  facet_grid(variable ~ ., scales = "free") + 
  scale_x_date()+ theme_minimal() + geom_smooth(method = "lm") +
  theme(legend.position = "none")
  
ggsave("compare_detrended_t_avg.png", p1, width = 15, height = 15, bg = "white")

# Fit a linear model for each group
models <- climate %>%
  group_by(variable) %>%
  do(model = lm(value ~ date, data = .))

# Extract the model summaries
model_summaries <- models %>%
  do(tidy(.$model)) %>%
  ungroup() %>%
  mutate(group = rep(models$variable, each = length(coef(models$model[[1]]))),
         p.value = format(round(p.value, 3), nsmall = 2))

# B) Example with precipitation

climate <- read.csv("Psum_SEZ.txt", sep = " ", header = TRUE)
climate <- select(climate, date, p_sum)
tmp_model <- lm(p_sum ~ seq(1:nrow(climate)), data = climate)
climate$p_sum_predicted <- predict(tmp_model)

climate$residual <- climate$p_sum - climate$p_sum_predicted
climate$detrended <- climate$residual / sd(climate$residual)
climate <- select(climate, date, p_sum, detrended)

climate <- melt(climate, id.vars = "date")
climate$date <- ymd(climate$date)


p1 <- ggplot(climate, aes(x = date, y = value,  colour = variable)) + geom_line() +
  facet_grid(variable ~ ., scales = "free") + 
  scale_x_date()+ theme_minimal() + geom_smooth(method = "lm") +
  theme(legend.position = "none")

ggsave("compare_detrended_p_sum.png", p1, width = 15, height = 15, bg = "white")

# Fit a linear model for each group
models <- climate %>%
  group_by(variable) %>%
  do(model = lm(value ~ date, data = .))

# Extract the model summaries
model_summaries <- models %>%
  do(tidy(.$model)) %>%
  ungroup() %>%
  mutate(group = rep(models$variable, each = length(coef(models$model[[1]]))),
         p.value = format(round(p.value, 3), nsmall = 2))



################################################################################
# Compare detrended and non-detrended results

chronology <- read.rwl("Sez.rwl")
detrended <- detrend(chronology, method = "Spline", nyrs = 33)
site_chron <- chron(detrended, prewhiten = TRUE)
residual <- site_chron[,2, drop = FALSE]

precipitation <- read.csv("Psum_SEZ.txt", sep = " ", header = TRUE)
precipitation <- select(precipitation, date, p_sum)
daily_precipitation <- data_transform(precipitation, format = "daily")

ds_daily <- daily_response(response = residual, env_data = daily_precipitation,
                           row_names_subset = TRUE, remove_insignificant = FALSE,
                           alpha = 0.05, day_interval = c(-150, 300),
                           lower_limit = 21, upper_limit = 180)

melted <- melt(data.frame(ds_daily$calculations))
melted$season_length <- seq(21, 180)
melted$variable <- as.numeric(gsub("X", "", melted$variable))
melted <- dplyr::filter(melted, !is.na(value))
melted$type <- "non-detrended"

#
  
ds_daily <- daily_response(response = residual, env_data = daily_precipitation,
                           row_names_subset = TRUE, remove_insignificant = FALSE,
                           dc_method = "SLD",
                           alpha = 0.05, day_interval = c(-150, 300),
                           lower_limit = 21, upper_limit = 180)

melted_dc <- melt(data.frame(ds_daily$calculations))
melted_dc$season_length <- seq(21, 180)
melted_dc$variable <- as.numeric(gsub("X", "", melted_dc$variable))
melted_dc <- dplyr::filter(melted_dc, !is.na(value))
melted_dc$type <- "detrended"

# Combine data frame for detrended and non-detrended data
binded <- rbind(melted_dc, melted)

p1 <- ggplot(binded,aes_(x = ~as.numeric(variable), y = ~season_length, fill = ~value)) +
  geom_tile() +
  xlab("Month") +
  facet_grid(type~.) +
  xlab("Day of Year") +
  ylab("Season Length") +
  scale_x_continuous(expand=c(0,0),breaks = c(150, 200, 250, 300, 350, 416, 466, 516, 566, 616), 
                                  labels = c("150*","200*", "250*", "300*", "350*", "50", "100", "150", "200", "250")) +
  scale_y_continuous(expand=c(0,0), breaks = seq(30, 270, by = 30)) +
  geom_vline(xintercept = 366, linetype="dashed", alpha = 0.5) +
  scale_fill_gradient2(low = "red", mid = "white", high = "blue", na.value = 'gray95', midpoint = 0, limits = c(-0.70, 0.70)) +
  theme_minimal() +  theme(axis.text = element_text(size = 12),
                           axis.title.y = element_text(size = 12), text = element_text(size = 12),
                           axis.title.x = element_blank(),
                           plot.title = element_text(size = 12),
                           panel.grid = element_blank(),
                           legend.title = element_blank(),
                           legend.position = "bottom", legend.key.width = unit(3, "line"),
                           panel.background = element_rect(fill = "gray95",
                                                           colour = "gray95",
                                                           size = 0.5, linetype = "solid"))

ggsave("Compare_detrended_non_detrended.png", p1, width = 9, height = 12, bg = "white")
