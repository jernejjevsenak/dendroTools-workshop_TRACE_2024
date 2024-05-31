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


# select window size / season_length (in days)
season_length <- 60

# select the number of years to be considered in each step
n_years <- 30

# define the first year
first_year <- 1950
last_year <- 2020


# we create a list to save output matrix from each step



for (x in c(40, 50, 60)){

  list_output <- list()
  place_holder <- 1
  
for (j in seq(1950, 1990)) {

  season_length <- x 
  
ds_daily_temp   <- daily_response(response = residual, env_data = daily_precipitation,
                                 row_names_subset = TRUE, remove_insignificant = FALSE,
                                 fixed_width = season_length,
                                 subset_years = c(j, j + n_years))

df_temp <- data.frame(ds_daily_temp$calculations)
df_temp$first_year <- j
df_temp$climate <- "Precipitation"

list_output[[place_holder]] <- df_temp
place_holder <- place_holder + 1

}




binded <- do.call(rbind, list_output)
binded <- melt(binded, id.vars = c("first_year", "climate"), variable.name = "season_start")

head(binded)

binded$season_start <- as.numeric(binded$season_start)
binded$first_year_label <- paste0(binded$first_year, "-",(binded$first_year + 30))

head(binded)

ggplot(binded, aes(x = season_start, y = first_year_label, fill = value)) + geom_tile() +
  scale_y_discrete(breaks = c(paste0(seq(1950, 1990, by = 5), "-", (seq(1950, 1990, by = 5)+30)))) +
  scale_x_continuous(expand = c(0,0), limits = c(1,322)) +
  scale_fill_gradient2() +
  xlab("Day of Year") + ylab("Period") + labs(fill = NULL)

ggsave(paste0("temporal_stability_daily_season",x, ".png"), width = 8, height = 8)

}



# Use critical threshold to remove insignificant correlations
critical_r <- function(n, alpha = .05) {
  df <- n - 2
  critical_t <- qt(alpha / 2, df, lower.tail = FALSE)
  critical_r_value <- sqrt((critical_t ^ 2) / ((critical_t ^ 2) + df))
  return(critical_r_value)
}

critical_r(40, 0.0001)

binded$value <- ifelse(abs(binded$value) < critical_r(n = 31, alpha = 0.05), NA, binded$value)
binded <- dplyr::filter(binded, !is.na(value))

ggplot(binded, aes(x = season_start, y = first_year_label, fill = value)) + geom_tile() +
  scale_y_discrete(breaks = c(paste0(seq(1950, 1990, by = 5), "-", (seq(1950, 1990, by = 5)+30)))) +
  scale_x_continuous(expand = c(0,0),limits = c(1,322)) +
  scale_fill_gradient2() +
  xlab("Day of Year") + ylab("Period") + labs(fill = NULL)

ggsave("temporal_stability_daily_filter.png", width = 8, height = 8)



#####################################################
#### Temporal stability for monthly correlations ####
#####################################################

# select window size / season_length (in months)
season_length <- 2

# select the number of years to be considered in each step
n_years <- 30

# define the first year
first_year <- 1950
last_year <- 2020


# we create a list to save output matrix from each step
list_output <- list()
place_holder <- 1

for (j in seq(1950, 1990)) {
  
  ds_monthly_temp   <- monthly_response(response = residual, env_data = monthly_precipitation,
                                    row_names_subset = TRUE, remove_insignificant = FALSE,
                                    fixed_width = season_length,
                                    subset_years = c(j, (j + n_years - 1)))
  
  df_temp <- data.frame(ds_monthly_temp$calculations)
  df_temp$first_year <- j
  df_temp$climate <- "Precipitation"
  
  list_output[[place_holder]] <- df_temp
  place_holder <- place_holder + 1
  
}


binded <- do.call(rbind, list_output)
binded <- melt(binded, id.vars = c("first_year", "climate"), variable.name = "season_start")
binded$season_start <- as.numeric(binded$season_start)
binded$first_year_label <- paste0(binded$first_year, "-",(binded$first_year + 30))

ggplot(binded, aes(x = season_start, y = first_year_label, fill = value)) + geom_tile() +
  scale_y_discrete(breaks = c(paste0(seq(1950, 1990, by = 5), "-", (seq(1950, 1990, by = 5)+30)))) +
  scale_x_continuous(expand = c(0,0),
                     breaks = seq(1,11),
                     labels = c("J","F","M","A","M","J","J","A","S","O","N")) +
  scale_fill_gradient2() +
  xlab("Day of Year") + ylab("Period") + labs(fill = NULL)

ggsave("temporal_stability_monthly.png", width = 8, height = 8)

# Use critical threshold to remove insignificant correlations
critical_r <- function(n, alpha = .05) {
  df <- n - 2
  critical_t <- qt(alpha / 2, df, lower.tail = FALSE)
  critical_r_value <- sqrt((critical_t ^ 2) / ((critical_t ^ 2) + df))
  return(critical_r_value)
}



binded$value <- ifelse(abs(binded$value) < critical_r(n = 31, alpha = 0.05), NA, binded$value)
binded <- dplyr::filter(binded, !is.na(value))

ggplot(binded, aes(x = season_start, y = first_year_label, fill = value)) + geom_tile() +
  scale_y_discrete(breaks = c(paste0(seq(1950, 1990, by = 5), "-", (seq(1950, 1990, by = 5)+30)))) +
  scale_x_continuous(expand = c(0,0),
                     limits = c(1,11), 
                     breaks = seq(1,11),
                     labels = c("J","F","M","A","M","J","J","A","S","O","N")) +
  scale_fill_gradient2() +
  xlab("Day of Year") + ylab("Period") + labs(fill = NULL)

ggsave("temporal_stability_monthly_filter.png", width = 8, height = 8)



