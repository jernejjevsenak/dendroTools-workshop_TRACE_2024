library(dplyr)
library(ggplot2)
library(rTG)
library(MLmetrics)

# Open data from rTG
data(data_trees)
data(parameters)

# run simulation 1 
simulation_1 <- XPSgrowth(data_trees, parameters,
                          ID_vars = c("Species", "Tissue", "Site", "Year", "Tree"),
                          fitting_method = c("brnn", "gam", "gompertz"),
                          search_initial_gom = TRUE,
                          fitted_save = TRUE, add_zeros = TRUE,
                          add_zeros_before = 75, post_process = TRUE)

plot(simulation_1)
summary(simulation_1)

# save simulation_1 output in separate data frame and remove manually added zeros
df_A <- simulation_1[[1]]
df_A <- mutate(df_A, width = ifelse(note == "added zero", NA, width))


# Improve labels for plotting purposes
df_A <- mutate(df_A, method = ifelse(method == "brnn", "BRNN", 
                                     ifelse(method == "gam", "GAM", 
                                            ifelse(method == "gompertz", "Gompertz", NA))))

# Sort conifer and broadlevaes
df_A$Species <- factor(df_A$Species, levels = c("FASY", "QUPU", "PCAB"))

# Create and save ggplot object
ggplot(df_A, aes(x = doy, y = width_pred, col = factor(Tree), group = key)) + 
  geom_line() + geom_point(df_A, mapping = aes(x = doy, y = width),alpha = 0.2) +
  facet_grid(Species + Tissue ~ method, scales = "free") +
  scale_colour_brewer(palette = "Dark2") + theme_minimal() +
  theme(legend.position = "none", text = element_text(size=18)) +
  ylab("Width [µm] / Cell Number") + xlab("Day of Year") +
  theme(panel.background = element_rect(fill = "gray97", colour = "gray80",
                                        size = 0.5, linetype = "solid")) 

ggsave("rTG_fit.png", width = 10, height = 10)

# Calculate RMSE, %RMSE and R2
group_by(dplyr::filter(df_A, !is.na(width)), Tissue, method, Species) %>%
                 summarise(RMSE = RMSE(width_pred, width),
                           rRMSE = RMSE/mean(width)*100,
                           R2_Score = MLmetrics::R2_Score(width, width_pred)) %>%
                 arrange(-RMSE) %>% arrange(Tissue, Species, method)



# Open data from rTG
data(data_trees)
data(parameters)

# Select one tree for XYLEM
data_trees <- dplyr::filter(data_trees, Species == "QUPU", Tissue == "XYLEM", Tree == 2)

# Simulation with no correction
simulation_A <- XPSgrowth(data_trees, parameters,
                          ID_vars = c("Species", "Tissue", "Site", "Year", "Tree"),
                          fitting_method = c("brnn", "gam"),
                          search_initial_gom = TRUE, add_zeros_before = 75,
                          fitted_save = FALSE, unified_parameters = FALSE,
                          add_zeros = FALSE, post_process = FALSE)

# Simulation with added zeros
simulation_B <- XPSgrowth(data_trees, parameters,
                          ID_vars = c("Species", "Tissue", "Site", "Year", "Tree"),
                          fitting_method = c("brnn", "gam"),
                          search_initial_gom = TRUE, add_zeros_before = 75,
                          fitted_save = FALSE, unified_parameters = FALSE,
                          add_zeros = TRUE, post_process = FALSE)

# Simulation with added zeros and post-correction
simulation_C <- XPSgrowth(data_trees, parameters,
                          ID_vars = c("Species", "Tissue", "Site", "Year", "Tree"),
                          fitting_method = c("brnn", "gam"),
                          search_initial_gom = TRUE, add_zeros_before = 75,
                          fitted_save = FALSE, unified_parameters = FALSE,
                          add_zeros = TRUE, post_process = TRUE)

# Extract the first elements from each simulation and give description
df_a <- simulation_A[[1]]; df_a$app <- "No correction" 
df_b <- simulation_B[[1]]; df_b$app <- "Added zeros"
df_c <- simulation_C[[1]]; df_c$app <- "Added zeros and\npost-correction"

# The next one is only to produce labels
df_abc <- simulation_B[[1]]; df_abc$method <- ifelse(df_abc$method == "gam", "GAM", "BRNN")

# Rbind a, b and c data frames and improve labels
temp_data <- rbind( df_a, df_b, df_c)
temp_data$app <- factor(temp_data$app, levels = c("No correction", "Added zeros", "Added zeros and\npost-correction"))
temp_data$method <- ifelse(temp_data$method == "gam", "GAM", "BRNN")

# create a data frame with labels 
label_df <- data.frame(ID = c("(a)","(b)","(c)","(d)","(e)","(f)"), 
                       method = c("BRNN", "BRNN","BRNN", "GAM","GAM", "GAM"), 
                       app = c("No correction", "Added zeros", "Added zeros and\npost-correction",
                               "No correction", "Added zeros", "Added zeros and\npost-correction"),
                       doy = 1, width_pred = 490)

label_df$app <- factor(label_df$app, levels = c("No correction", "Added zeros", "Added zeros and\npost-correction"))

# Create ggplot object
ggplot(temp_data, aes(x = doy, y = width_pred)) +
  geom_line() + facet_grid(app~method) + theme_bw() +
  geom_point(df_abc, mapping = aes(x = doy, y = width, alpha = note)) +
  scale_alpha_discrete(range=c(0.1, 1)) + guides(alpha = "none") +
  ylab("Xylem Width [µm]") + xlab("Day of Year") + theme_minimal() +
  theme(panel.background = element_rect(fill = "gray97", colour = "gray80", 
                                        size = 0.5, linetype = "solid"),
                                        text = element_text(size = 14)) +
  geom_text(data = label_df, mapping = aes(label = ID))


















# open dendrometer data
data("data_dendrometers")

# apply the XPSgrowth() with the post-processing algorithm
sim1 <- XPSgrowth(data_dendrometers, unified_parameters = TRUE,
                  ID_vars = c("site", "species", "year", "tree"),
                  fitting_method = c("brnn", "gam", "gompertz"), 
                  brnn_neurons = 2, gam_k = 9, gam_sp = 0.5, 
                  search_initial_gom = TRUE,
                  add_zeros = FALSE, post_process = TRUE)

# extract the fitted values and add label "With post-process"
simulation_A <- sim1$fitted
simulation_A$label <- "With post-process"

# apply the XPSgrowth() without the post-processing algorithm
sim2 <- XPSgrowth(data_dendrometers, unified_parameters = TRUE,
                  ID_vars = c("site", "species", "year", "tree"),
                  fitting_method = c("brnn", "gam", "gompertz"), 
                  brnn_neurons = 2, gam_k = 9, gam_sp = 0.5, 
                  search_initial_gom = TRUE, add_zeros = FALSE, 
                  post_process = FALSE)

# extract the fitted values and add label "Without post-process"
simulation_B <- sim2$fitted
simulation_B$label <- "Without post-process"

# rbind both simulations and improve method labels
combined <- rbind(simulation_A, simulation_B)
combined <- dplyr::mutate(combined, 
                          method = ifelse(method == "brnn", "BRNN", 
                                         ifelse(method == "gam", "GAM", 
                                               ifelse(method == "gompertz", "Gompertz", NA))))
# create ggplot
ggplot(combined, aes(x = doy, y = width)) + geom_point() +
  geom_line(aes(x = doy, y = width_pred), col = "red") +
  facet_grid(species + label ~ method) + scale_colour_brewer(palette = "Dark2") + theme_minimal() +
  theme(legend.position = "none", text = element_text(size=18)) +
  ylab("Width") + xlab("Day of Year") +
  theme(panel.background = element_rect(fill = "gray97",
                                        colour = "gray80",
                                        size = 0.5, linetype = "solid"))

ggsave("dendrometer_example.png", height = 9, width = 10)
