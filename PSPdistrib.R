library(tidyverse)
library(fitdistrplus)
library(foreach)
library(doParallel)

# ==========================================
# 1. LOAD AND PRE-PROCESS DATA
# ==========================================
raw_data <- read_csv("forestry_data.csv")

tree_data <- raw_data %>%
  mutate(
    TreeNo = ifelse(is.na(TreeNo), TREE_KEY, TreeNo),
    Fork = ifelse(is.na(Fork) | Fork == "" | Fork == "A", "A", Fork),
    Age = as.numeric(difftime(as.Date(MeasDate), as.Date(Planted), units = "days")) / 365.25
  ) %>%
  rename(
    Region       = Code,
    Plot         = PLOT_KEY,
    PlantedDate  = Planted,
    MeasureDate  = MeasDate
  ) %>%
  dplyr::select(
    Species, Region, Plot, PlantedDate, MeasureDate, TreeNo, Fork, DBH, TotalHt,
    PlotArea, Plot_Id, MEAS_KEY, THIN_KEY, PRUNE_KEY, ExpNo, SubExp, PlotNo,
    SubPlot, Forest, Compartment, Stand, RowSpacing, TreeSpacing, Aspect, SiteIndex, Age
  )

# Filter for missing values and restrict strictly to P. radiata
tree_data <- tree_data %>%
  filter(
    !is.na(Region),
    !is.na(DBH),
    !is.na(TotalHt),
    !is.na(Age),
    Species == "P.RAD"
  )

# ==========================================
# 2. STANDARDIZE REGIONS AND SCALE (0 TO 1)
# ==========================================
processed_data <- tree_data %>%
  mutate(
    Decade  = paste0(floor(as.numeric(format(as.Date(MeasureDate), "%Y")) / 10) * 10, "s"),
    Age_Bin = round(Age)
  ) %>%
  # Map raw codes strictly to your official NZ names list
  mutate(Region = case_when(
    Region == "AK" ~ "Auckland",
    Region == "BP" ~ "Bay of Plenty",
    Region == "CY" ~ "Canterbury",
    Region == "WD" ~ "West Coast",
    Region == "FR" ~ "West Coast",  
    Region == "JN" ~ "Gisborne",    
    Region == "NN" ~ "Nelson",
    Region == "BC" ~ "Marlborough", 
    Region == "SD" ~ "Southland",
    Region == "SR" ~ "Southland",   
    Region == "ST" ~ "SN Island",   
    Region == "PR" ~ "Waikato",     
    Region == "RO" ~ "Taupo",       
    TRUE           ~ NA_character_  
  )) %>%
  # Filter out ages > 50 and invalid regions
  filter(
    Age_Bin <= 50,
    !is.na(Region)
  ) %>%
  # Calculate cohort-specific maximums to scale data between 0 and 1
  group_by(Region, Age_Bin) %>%
  mutate(
    Max_DBH_Cohort = max(DBH, na.rm = TRUE),
    Max_HT_Cohort  = max(TotalHt, na.rm = TRUE),
    DBH_Rel        = DBH / Max_DBH_Cohort,
    HT_Rel         = TotalHt / Max_HT_Cohort
  ) %>%
  ungroup()

rm(raw_data)
rm(tree_data)

# ==========================================
# 3. COMPETITIVE DISTRIBUTION SELECTION (AIC)
# ==========================================
get_best_dist <- function(data_vector) {
  clean_data <- data_vector[!is.na(data_vector) & data_vector > 0]
  
  if (length(clean_data) < 5) {
    return(tibble(best_dist = "Insufficient Data", w_aic = NA, g_aic = NA, ln_aic = NA))
  }
  
  fit_w  = try(suppressWarnings(fitdist(clean_data, "weibull", method = "mle")), silent = TRUE)
  fit_g  = try(suppressWarnings(fitdist(clean_data, "gamma", method = "mle")), silent = TRUE)
  fit_ln = try(suppressWarnings(fitdist(clean_data, "lnorm", method = "mle")), silent = TRUE)
  
  aic_w  <- if (inherits(fit_w, "fitdist")) fit_w$aic else Inf
  aic_g  <- if (inherits(fit_g, "fitdist")) fit_g$aic else Inf
  aic_ln <- if (inherits(fit_ln, "fitdist")) fit_ln$aic else Inf
  
  if (all(c(aic_w, aic_g, aic_ln) == Inf)) {
    return(tibble(best_dist = "Fit Failed", w_aic = NA, g_aic = NA, ln_aic = NA))
  }
  
  aics <- c(Weibull = aic_w, Gamma = aic_g, Lognormal = aic_ln)
  best <- names(aics)[which.min(aics)]
  
  return(tibble(
    best_dist = best, 
    w_aic     = ifelse(aic_w == Inf, NA, aic_w), 
    g_aic     = ifelse(aic_g == Inf, NA, aic_g), 
    ln_aic    = ifelse(aic_ln == Inf, NA, aic_ln)
  ))
}

# Run the AIC analysis sequentially (very quick compared to plotting)
dbh_distribution_summary <- processed_data %>%
  group_by(Region, Decade, Age_Bin) %>%
  do(get_best_dist(.$DBH_Rel)) %>%
  ungroup() %>%
  rename_with(~ paste0("dbh_", .), c(best_dist, w_aic, g_aic, ln_aic))

ht_distribution_summary <- processed_data %>%
  group_by(Region, Decade, Age_Bin) %>%
  do(get_best_dist(.$HT_Rel)) %>%
  ungroup() %>%
  rename_with(~ paste0("ht_", .), c(best_dist, w_aic, g_aic, ln_aic))

final_distribution_summary <- full_join(
  dbh_distribution_summary, 
  ht_distribution_summary, 
  by = c("Region", "Decade", "Age_Bin")
)

write_csv(final_distribution_summary, "best_fit_distributions_summary.csv")

# ==========================================
# 4. PREPARE THE MODEL PARAMETERS AND CURVES
# ==========================================
weibull_models_ht <- processed_data %>%
  filter(!is.na(HT_Rel) & HT_Rel > 0) %>%
  group_by(Region, Decade, Age_Bin) %>%
  filter(n() >= 5) %>% 
  summarise(
    fit_obj = list(tryCatch(
      suppressWarnings(fitdist(HT_Rel, "weibull", method = "mle")),
      error = function(e) NULL
    )),
    .groups = "drop"
  ) %>%
  filter(!map_lgl(fit_obj, is.null)) %>%
  mutate(
    shape = map_dbl(fit_obj, ~ .$estimate["shape"]),
    scale = map_dbl(fit_obj, ~ .$estimate["scale"])
  )

curve_points_ht <- weibull_models_ht %>%
  rowwise(Region, Decade, Age_Bin) %>%
  reframe(
    HT_Rel  = seq(0, 1, length.out = 100),
    Density = dweibull(HT_Rel, shape = shape, scale = scale)
  )

weibull_models_dbh <- processed_data %>%
  filter(!is.na(DBH_Rel) & DBH_Rel > 0) %>%
  group_by(Region, Decade, Age_Bin) %>%
  filter(n() >= 5) %>% 
  summarise(
    fit_obj = list(tryCatch(
      suppressWarnings(fitdist(DBH_Rel, "weibull", method = "mle")),
      error = function(e) NULL
    )),
    .groups = "drop"
  ) %>%
  filter(!map_lgl(fit_obj, is.null)) %>%
  mutate(
    shape = map_dbl(fit_obj, ~ .$estimate["shape"]),
    scale = map_dbl(fit_obj, ~ .$estimate["scale"])
  )

curve_points_dbh <- weibull_models_dbh %>%
  rowwise(Region, Decade, Age_Bin) %>%
  reframe(
    DBH_Rel = seq(0, 1, length.out = 100),
    Density = dweibull(DBH_Rel, shape = shape, scale = scale)
  )

# ==========================================
# 5. INITIALIZE PARALLEL CLUSTER FOR PLOTS
# ==========================================
cores <- parallel::detectCores() - 1
cl    <- makeCluster(cores)
registerDoParallel(cl)

# FIXED: Export the required curve coordinates alongside the raw points
clusterEvalQ(cl, {
  library(tidyverse)
})
clusterExport(cl, c("processed_data", "curve_points_ht", "curve_points_dbh"))

# ==========================================
# 6. PARALLEL HEIGHT PLOTTING ENGINE
# ==========================================
unique_ages_ht <- unique(weibull_models_ht$Age_Bin)

foreach(a = unique_ages_ht, .packages = c("tidyverse")) %dopar% {
  dir_path <- paste0("Age_", a)
  if (!dir.exists(dir_path)) dir.create(dir_path)
  
  p_sub_age <- processed_data %>% filter(Age_Bin == a)
  c_sub_age <- curve_points_ht %>% filter(Age_Bin == a)
  regions_in_age <- unique(c_sub_age$Region)
  
  for (r in regions_in_age) {
    p_sub <- p_sub_age %>% filter(Region == r)
    c_sub <- c_sub_age %>% filter(Region == r)
    
    p <- ggplot() +
      geom_histogram(data = p_sub, aes(x = HT_Rel, y = after_stat(density), fill = Decade), 
                     binwidth = 0.02, position = "identity", alpha = 0.3) +
      geom_line(data = c_sub, aes(x = HT_Rel, y = Density, color = Decade), size = 1) +
      coord_cartesian(xlim = c(0, 1)) + 
      scale_x_continuous(expand = c(0, 0)) +
      labs(title = paste("Region:", r, "- Age:", a, "- Relative Height Profile"),
           x = "Relative Tree Height (Proportion of Max)", y = "Probability Density") +
      theme_minimal()
    
    ggsave(file.path(dir_path, paste0(r, "_relative_height.png")), plot = p, width = 8, height = 5, dpi = 300)
  }
}

# ==========================================
# 7. PARALLEL DBH PLOTTING ENGINE
# ==========================================
unique_ages_dbh <- unique(weibull_models_dbh$Age_Bin)

foreach(a = unique_ages_dbh, .packages = c("tidyverse")) %dopar% {
  dir_path <- paste0("Age_", a)
  if (!dir.exists(dir_path)) dir.create(dir_path)
  
  p_sub_age <- processed_data %>% filter(Age_Bin == a)
  c_sub_age <- curve_points_dbh %>% filter(Age_Bin == a)
  regions_in_age <- unique(c_sub_age$Region)
  
  for (r in regions_in_age) {
    p_sub <- p_sub_age %>% filter(Region == r)
    c_sub <- c_sub_age %>% filter(Region == r)
    
    p <- ggplot() +
      geom_histogram(data = p_sub, aes(x = DBH_Rel, y = after_stat(density), fill = Decade), 
                     binwidth = 0.02, position = "identity", alpha = 0.3) +
      geom_line(data = c_sub, aes(x = DBH_Rel, y = Density, color = Decade), size = 1) +
      coord_cartesian(xlim = c(0, 1)) +
      scale_x_continuous(expand = c(0, 0)) +
      labs(title = paste("Region:", r, "- Age:", a, "- Relative DBH Profile"),
           x = "Relative DBH (Proportion of Max)", y = "Probability Density") +
      theme_minimal()
    
    ggsave(file.path(dir_path, paste0(r, "_relative_dbh.png")), plot = p, width = 8, height = 5, dpi = 300)
  }
}

# ==========================================
# 8. SHUTDOWN PARALLEL ENGINE & EXPORT TRENDS
# ==========================================
stopCluster(cl)

# Save the metric statistics to file
weibull_param_summary <- full_join(
  weibull_models_dbh %>% dplyr::select(Region, Decade, Age_Bin, dbh_shape = shape, dbh_scale = scale),
  weibull_models_ht  %>% dplyr::select(Region, Decade, Age_Bin, ht_shape  = shape, ht_scale  = scale),
  by = c("Region", "Decade", "Age_Bin")
)

write_csv(weibull_param_summary, "weibull_parameters_trend_analysis.csv")

