library("dplR")
library("dendroTools")
library("lubridate")
library("ggplot2")
library("dplyr")
library("reshape2")

# Open TRWi data
temp_rwl <- read.crn("my_chron.crn")
temp_rwl$samp.depth <- NULL # remove the sample depth if present

# Open climate data
Prec <- read.table("my_chron_Psum.csv", sep = ",", header = TRUE)

Prec$date <- ymd(paste0(Prec$Y,"-", Prec$M, "-", Prec$D))

Prec <- dplyr::mutate(Prec, year = year(date), doy = yday(date)) %>%
  select(year, doy, Prec)

Prec <- reshape2::dcast(data = Prec, formula = year~doy)
Prec <- years_to_rownames(Prec, "year")

ds_p <- daily_response(response = temp_rwl, env_data = Prec, method = "cor",
               row_names_subset = TRUE, previous_year = FALSE, lower_limit = 30,
               upper_limit = 90, reference_window = "end")

output <- data.frame(ds_p$calculations)
output <- melt(output)
output$season_length <- seq(30, 90)

ggplot(output,aes_(x = ~as.numeric(variable), y = ~season_length, fill = ~value)) +
  geom_tile() +
  xlab("Day of Year") +
  ylab("Season Length") +
  scale_x_continuous(expand=c(0,0),breaks = c(20, 150, 300, 450, 584, 725), labels = c("jan*", "may*", "oct*", "MAR", "AUG", "NOV")) +
  scale_y_continuous(expand=c(0,0), breaks = seq(30, 270, by = 30)) +
  geom_vline(xintercept = 366) +
  scale_fill_gradient2(low = "red", mid = "white", high = "blue", na.value = 'gray97', midpoint = 0, limits = c(-0.8, 0.8)) +
  theme_minimal() +  theme(axis.text = element_text(size = 10),
                           axis.title.y = element_text(size = 18), text = element_text(size = 18),
                           axis.title.x = element_blank(),
                           plot.title = element_text(size = 16),
                           legend.title = element_blank(),
                           legend.position = "bottom", legend.key.width = unit(3, "line"),
                           panel.background = element_rect(fill = "gray97",
                                                           colour = "gray80",
                                                           size = 0.5, linetype = "solid"))

ggsave("Precipitation_example_daily.png", height = 7, width = 10)
