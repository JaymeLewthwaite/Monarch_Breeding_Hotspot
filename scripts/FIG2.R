############################################################################
###   SCRIPT FOR ANALYSIS AND FIGURE 2 in LEWTHWAITE ET AL. 2026    ########
###   ANALYSIS CODED BY FR, REVIEWED AND ENHANCED BY JL		    ########
############################################################################

# libraries 
library(glmmTMB)
library(lme4)
library(dplyr)
library(ggplot2)
library(ggeffects)
library(lubridate)
library(DHARMa)
library(performance)
library(maps)
library(patchwork)
library(tidyverse)
library(sf)
library(purrr)
library(tibble)
library(bbmle)   

############################################################################
###   DATA PREP    #########################################################
############################################################################

data <- read.csv("breeding_df.csv")

# re-organize to have the grid cells at the beginning of the dataset.
# remove info on locality, country etc
data <- data[, c(1, 49, 2, 50, 51, 3:48)]
data <- data[, -c(9, 10, 11)]

# remove locations with 0 counts for the stems, as it breaks the offset term
data <- subset(
  data,
  !is.na(LarvaeCount) &
    !is.na(StemCount) &
    StemCount > 0 &
    !is.na(YearCollected) &
    !is.na(Ecoregion))

## open map for crs
grid_50 <- st_read("grid_50km/grid_50km.shp")

# random effect variables must be coded as factors
data$YearCollected <- as.factor(data$YearCollected)
data$Ecoregion     <- as.factor(data$Ecoregion)

# remove old columns
data <- data[-c(2:5)]

data <- data %>%
  mutate(ObservationDate = as.Date(ObservationDate, format = "%d/%m/%Y"),
    YearCollected   = as.factor(YearCollected),
    Ecoregion       = as.factor(Ecoregion),
    WeekOfYear      = as.factor(isoweek(ObservationDate))
  ) %>%
  filter(StemCount > 0)

# multiple locations with the same coordinates and in the same year
set.seed(123)   

# select for each location with multiple observations in the same year, 
# only one observation

# maximum value, for sensitivity analysis
data_max <- data %>%
  group_by(DecimalLongitude, DecimalLatitude, YearCollected) %>%
  slice_max(LarvaeCount, n = 1, with_ties = FALSE) %>%
  ungroup()

# random value, used in main analysis
data <- data %>%
  group_by(DecimalLongitude, DecimalLatitude, YearCollected) %>%
  slice_sample(n = 1) %>%
  ungroup()

# points sf object in WGS84
# project to Albers
points <- st_as_sf(data,coords = c("DecimalLongitude", "DecimalLatitude"),crs = 4326)
points <- st_transform(points, st_crs(grid_50))

############################################################################
###   SET UP GRID    #######################################################
############################################################################

library(h3jsr)

# reproject to WGS84 (required by h3jsr)
points_4326 <- st_transform(points, crs = 4326)

# assign finest resolution directly from point coordinates
points_4326$h3_res7 <- point_to_cell(points_4326, res = 7)  # ~ 26 km2

# derive all coarser resolutions from res7 to guarantee perfect nesting
points_4326$h3_res6 <- h3jsr::get_parent(points_4326$h3_res7, res = 6)  # ~ 91 km2
points_4326$h3_res5 <- h3jsr::get_parent(points_4326$h3_res7, res = 5)  # ~ 252 km2
points_4326$h3_res4 <- h3jsr::get_parent(points_4326$h3_res7, res = 4)  # ~ 737 km2
points_4326$h3_res3 <- h3jsr::get_parent(points_4326$h3_res7, res = 3)  # ~ 737 km2

# # verify nesting 
# cat("res7 nested in res6:", all(h3jsr::get_parent(points_4326$h3_res7, res = 6) == points_4326$h3_res6), "\n")
# cat("res6 nested in res5:", all(h3jsr::get_parent(points_4326$h3_res6, res = 5) == points_4326$h3_res5), "\n")
# cat("res5 nested in res4:", all(h3jsr::get_parent(points_4326$h3_res5, res = 4) == points_4326$h3_res4), "\n")
# cat("res4 nested in res3:", all(h3jsr::get_parent(points_4326$h3_res4, res = 3) == points_4326$h3_res3), "\n")

# copy IDs back to projected points
points$h3_res7 <- points_4326$h3_res7
points$h3_res6 <- points_4326$h3_res6
points$h3_res5 <- points_4326$h3_res5
points$h3_res4 <- points_4326$h3_res4
points$h3_res3 <- points_4326$h3_res3

# check no NAs 
# colSums(is.na(st_drop_geometry(points)[, c("h3_res7", "h3_res6", "h3_res5", "h3_res4")]))

# hex fishnet for plots
hex_res7 <- cell_to_polygon(unique(points_4326$h3_res7), simple = FALSE)
hex_res6 <- cell_to_polygon(unique(points_4326$h3_res6), simple = FALSE)
hex_res5 <- cell_to_polygon(unique(points_4326$h3_res5), simple = FALSE)
hex_res4 <- cell_to_polygon(unique(points_4326$h3_res4), simple = FALSE)
hex_res3 <- cell_to_polygon(unique(points_4326$h3_res3), simple = FALSE)

# reproject hex polygons to match points CRS
hex_res7 <- st_transform(hex_res7, st_crs(points))
hex_res6 <- st_transform(hex_res6, st_crs(points))
hex_res5 <- st_transform(hex_res5, st_crs(points))
hex_res4 <- st_transform(hex_res4, st_crs(points))
hex_res3 <- st_transform(hex_res3, st_crs(points))

# grid for supplementary material
hex_grid<-ggplot() +
  geom_sf(data = hex_res3, fill = NA, color = "darkred",    linewidth = 0.6) +
  geom_sf(data = hex_res4, fill = NA, color = "red",    linewidth = 0.6) +
  geom_sf(data = hex_res5, fill = NA, color = "orange", linewidth = 0.4) +
  geom_sf(data = hex_res6, fill = NA, color = "gold",   linewidth = 0.2) +
  geom_sf(data = points,   size = 0.5) +
  theme_minimal()

ggsave("diagnostics_and_checks/Plot - Map of nested, hexagonal spatial units.png",
       plot   = hex_grid,
       width  = 30,
       height = 20,
       dpi    = 300,
       bg     = "white")

# assign IDs to hex grid
hex_res7 <- hex_res7 %>% mutate(hex_id_res7 = h3_address)
hex_res6 <- hex_res6 %>% mutate(hex_id_res6 = h3_address)
hex_res5 <- hex_res5 %>% mutate(hex_id_res5 = h3_address)
hex_res4 <- hex_res4 %>% mutate(hex_id_res4 = h3_address)
hex_res3 <- hex_res3 %>% mutate(hex_id_res3 = h3_address)

# copy H3 indices already in points to named ID columns ---
points$hex_id_res7 <- points$h3_res7
points$hex_id_res6 <- points$h3_res6
points$hex_id_res5 <- points$h3_res5
points$hex_id_res4 <- points$h3_res4
points$hex_id_res3 <- points$h3_res3

# generate data for analysis
data_model <- st_drop_geometry(points[,c(1:4,13:16,44,50:54)])

# Extract latitude from the geometry before dropping it
points$latitude <- st_coordinates(points)[, 2]

data_model <- st_drop_geometry(points[, c(1:4, 13:16, 44, 50:54)]) |>
  bind_cols(latitude = points$latitude)

# replace larval counts and milkweed stems 
data_model_max <- data_model
data_model_max$LarvaeCount <- data_max$LarvaeCount
data_model_max$StemCount <- data_max$StemCount
data_model_max$JulianDay <- data_max$JulianDay #JL added this line
data_model_max$SurveyID <- factor(data_max$SurveyID) #JL added this line


############################################################################
###   SUPPL MAT   ##########################################################
############################################################################

###
### CHECK NUMBER OF POINTS PER HEXAGON
###
summary_res3 <- data_model %>%
  group_by(hex_id_res3) %>%
  summarise(
    n_points = n(),
    n_child = n_distinct(hex_id_res4),
    .groups = "drop"
  ) %>%
  mutate(resolution = "Res3")

summary_res4 <- data_model %>%
  group_by(hex_id_res4) %>%
  summarise(
    n_points = n(),
    n_child = n_distinct(hex_id_res5),
    .groups = "drop"
  ) %>%
  mutate(resolution = "Res4")

summary_res5 <- data_model %>%
  group_by(hex_id_res5) %>%
  summarise(
    n_points = n(),
    n_child = n_distinct(hex_id_res6),
    .groups = "drop"
  ) %>%
  mutate(resolution = "Res5")

summary_res6 <- data_model %>%
  group_by(hex_id_res6) %>%
  summarise(
    n_points = n(),
    n_child = n_distinct(hex_id_res7),
    .groups = "drop"
  ) %>%
  mutate(resolution = "Res6")

summary_res7 <- data_model %>%
  group_by(hex_id_res7) %>%
  summarise(
    n_points = n(),
    .groups = "drop"
  ) %>%
  mutate(
    n_child = NA,
    resolution = "Res7")

hex_summary <- bind_rows(
  summary_res3,
  summary_res4,
  summary_res5,
  summary_res6,
  summary_res7)

hex_summary$resolution <- factor(
  hex_summary$resolution,
  levels = c("Res3", "Res4", "Res5", "Res6", "Res7"),
  labels = c("12,500 km²","750 km²",
    "250 km²","90 km²", "25 km²"))

hex_summary_plot<-ggplot(hex_summary,
       aes(x = log2(n_points))) +
  geom_histogram(binwidth = 1,
                 fill = "steelblue",
                 colour = "black") +
  facet_wrap(~resolution,
             scales = "free_y",
             ncol = 3) +
  scale_x_continuous(
    name = "Number of points per hexagon",
    breaks = log2(c(1, 2, 4, 8, 16, 32, 64, 128, 256, 512)),
    labels = c(1, 2, 4, 8, 16, 32, 64, 128, 256, 512)
  ) +
  theme_bw(base_size = 13)

ggsave("diagnostics_and_checks/Plot - Number of points per hexagon.png",
       plot   = hex_summary,
       width  = 30,
       height = 20,
       dpi    = 300,
       bg     = "white")

# table for supplementary material
multi_obs_summary <- hex_summary %>%
  group_by(resolution) %>%
  summarise(
    `Total # of hexgrid cells` = n(),
    `Total # of hexgrid cells with > 1 observation` = sum(n_points > 1),
    `Proportion of hexgrid cells with > 1 observation` =
      paste0(round(100 * mean(n_points > 1)), "%"),
    .groups = "drop"
  ) %>%
  rename(Resolution = resolution) %>%
  arrange(desc(Resolution))

# View table
print(multi_obs_summary)

# Export table as CSV
write.csv(
  multi_obs_summary,
  "diagnostics_and_checks/Proportion_hex_grids_with_multiple_observations.csv",
  row.names = FALSE
)

multi_obs_summary

############################################################################
###   MODELS   #############################################################
############################################################################

###
### TEST OF ALTERNATIVE MODELS
###

# set up main model formula
fixed  <- LarvaeCount ~ 1 +
  poly(JulianDay, 2) * scale(latitude) +
  (1 | SurveyID)+
  (1 | hex_id_res3/hex_id_res4/hex_id_res5/hex_id_res6/hex_id_res7) +
  (1 | Ecoregion) +
  (1 | YearCollected)

offset_arg <- quote(log(StemCount))

# alternative distributional assumptions
m_pois <- glmmTMB(fixed, family = poisson(),   offset = log(StemCount), data = data_model)

m_nb   <- glmmTMB(fixed, family = nbinom2(),   offset = log(StemCount), data = data_model)

m_zip  <- glmmTMB(fixed, family = poisson(),   offset = log(StemCount),
                  ziformula = ~1, data = data_model)

m_zinb <- glmmTMB(fixed, family = nbinom2(),   offset = log(StemCount),
                  ziformula = ~1, data = data_model)

# check support based on AIC
model_comparison<-AICtab(m_pois, m_nb, m_zip, m_zinb,
       weights = TRUE, delta = TRUE, base = TRUE)
# best model: negative binomial with no zero inflation term

# Convert to a formatted table
model_comparison_table <- as.data.frame(model_comparison) %>%
  tibble::rownames_to_column("Model") %>%
  mutate(
    Model = recode(
      Model,
      "m_nb"   = "Negative binomial",
      "m_zinb" = "Negative binomial with zero inflation",
      "m_pois" = "Poisson",
      "m_zip"  = "Poisson with zero inflation"
    ),
    AIC = round(as.numeric(AIC), 1),
    `Delta AIC` = round(as.numeric(dAIC), 1),
    df = as.integer(df),
    AICw = case_when(
      is.na(weight) ~ NA_character_,
      as.numeric(weight) < 0.01 ~ "< 0.01",
      TRUE ~ sprintf("%.2f", as.numeric(weight))
    )
  ) %>%
  select(Model, AIC, `Delta AIC`, df, AICw) %>%
  arrange(AIC)

# View table
print(model_comparison_table)

# Export table
write.csv(model_comparison_table,
  "diagnostics_and_checks/AIC_model_comparison.csv",
  row.names = FALSE,
  na = "NA")

# check model diagnostics with DHARMa
sim_diagnostics <- function(model, model_name, n_sim = 500) {
  sim <- simulateResiduals(model, n = n_sim, plot = FALSE)
  
  par(mfrow = c(1, 3))
  plot(sim, main = model_name)
  plotResiduals(sim, data_model$JulianDay,  xlab = "Julian Day",  main = "vs Julian Day")
  plotResiduals(sim, data_model$latitude,   xlab = "Latitude",   main = "vs Latitude")
  par(mfrow = c(1, 1))
  
  list(ks_test     = testUniformity(sim),
    dispersion     = testDispersion(sim),
    zero_inflation = testZeroInflation(sim),
    outliers       = testOutliers(sim))
}

diag_pois  <- sim_diagnostics(m_pois,  "Poisson")
diag_nb    <- sim_diagnostics(m_nb,    "NB2 (baseline)")
diag_zip   <- sim_diagnostics(m_zip,   "ZIP")
diag_zinb  <- sim_diagnostics(m_zinb,  "ZINB")

extract_pvals <- function(diag, name) {
  data.frame(
    model          = name,
    KS_uniformity  = diag$ks_test$p.value,
    dispersion_p   = diag$dispersion$p.value,
    zero_infl_p    = diag$zero_inflation$p.value,
    outlier_p      = diag$outliers$p.value
  )
}

pval_table <- do.call(rbind, list(
  extract_pvals(diag_pois,  "Poisson"),
  extract_pvals(diag_nb,    "NB2"),
  extract_pvals(diag_zip,   "ZIP"),
  extract_pvals(diag_zinb,  "ZINB")
))

#print(pval_table, digits = 3)
## nb model selected for further analysis. 
## while there is a significant departures from the simulated residual distribution,
## according to DHARMa, visual inspection of QQ plot suggests modest deviations 
## from the expected distribution. Given the sample size, tests are detecting 
## small discrepancies that do not invalidate the overall model inferences

##
## FINAL MODEL
##

m_pheno_lat <- glmmTMB(
  LarvaeCount ~ 1 +
    poly(JulianDay, 2) * scale(latitude) +
    (1 | SurveyID)+
    (1 | hex_id_res3/hex_id_res4/hex_id_res5/hex_id_res6/hex_id_res7) +
    (1 | Ecoregion) +
    (1 | YearCollected), 
  family = nbinom2(),
  offset = log(StemCount),
  data   = data_model)

# sensitivity to selection of random observation in pseudoreplicated sites (replaced with maximum observation)
m_pheno_lat_max <- glmmTMB(
  LarvaeCount ~ 1 +
    poly(JulianDay, 2) * scale(latitude) +
    (1 | SurveyID)+
    (1 | hex_id_res3/hex_id_res4/hex_id_res5/hex_id_res6/hex_id_res7) +
    (1 | Ecoregion) +
    (1 | YearCollected), 
  family = nbinom2(),
  offset = log(StemCount),
  data   = data_model_max)

# check model
summary(m_pheno_lat)
summary(m_pheno_lat_max)

saveRDS(m_pheno_lat, "model_outputs/m_pheno_lat.rds")

# check whether the model outcome is sensitive to random selection of one 
# sample in years with multiple observations from the same location
compare_re <- function(mod1, mod2) {
  
  re1 <- ranef(mod1)$cond
  re2 <- ranef(mod2)$cond
  
  common <- intersect(names(re1), names(re2))
  
  out <- lapply(common, function(grp) {
    
    d1 <- data.frame(
      level = rownames(re1[[grp]]),
      re1 = re1[[grp]][, 1]
    )
    
    d2 <- data.frame(
      level = rownames(re2[[grp]]),
      re2 = re2[[grp]][, 1]
    )
    
    d <- merge(d1, d2, by = "level")
    
    data.frame(
      Group = grp,
      Levels = nrow(d),
      Correlation = round(cor(d$re1, d$re2), 2)
    )
  })
  
  do.call(rbind, out)
}

# Generate formatted table
re_comparison <- compare_re(m_pheno_lat, m_pheno_lat_max) %>%
  mutate(
    `Random effect` = recode(
      Group,
      "hex_id_res7:hex_id_res6:hex_id_res5:hex_id_res4:hex_id_res3" = "25 km2 grid",
      "hex_id_res6:hex_id_res5:hex_id_res4:hex_id_res3" = "90 km2 grid",
      "hex_id_res5:hex_id_res4:hex_id_res3" = "250 km2 grid",
      "hex_id_res4:hex_id_res3" = "750 km2 grid",
      "hex_id_res3" = "12,500 km2 grid",
      "SurveyID" = "Survey ID",
      "Ecoregion" = "Ecoregion",
      "YearCollected" = "Year of sampling"
    )
  ) %>%
  select(
    `Random effect`,
    `Number of levels` = Levels,
    Correlation
  ) %>%
  mutate(
    `Random effect` = factor(
      `Random effect`,
      levels = c(
        "25 km2 grid",
        "90 km2 grid",
        "250 km2 grid",
        "750 km2 grid",
        "12,500 km2 grid",
        "Survey ID",
        "Ecoregion",
        "Year of sampling"
      )
    )
  ) %>%
  arrange(`Random effect`) %>%
  mutate(`Random effect` = as.character(`Random effect`))

# View table
print(re_comparison)

# Export as CSV
write.csv(
  re_comparison,
  "diagnostics_and_checks/Random_effect_correlations.csv",
  row.names = FALSE
)

# check phenology 
# latitudes for prediction: mean ± 1 SD 
lat_mean <- mean(data_model$latitude, na.rm = TRUE)
lat_sd   <- sd(data_model$latitude, na.rm = TRUE)

pred <- ggpredict(m_pheno_lat,
  terms = c("JulianDay [all]",
    sprintf("latitude [%.2f, %.2f, %.2f]",
            lat_mean - lat_sd,
            lat_mean,
            lat_mean + lat_sd)))

ggplot(pred, aes(x = x, y = predicted,
           colour = group, fill = group)) +
  geom_line(size = 1.2) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
              alpha = 0.15, colour = NA) +
  labs(x = "Julian Day",
    y = "Predicted Larvae Count",
    colour = "Latitude",
    fill = "Latitude") +
  theme_bw()
# at higher latitude, number of caterpillars increased more at the later dates 
# consistent with our expectation that sites in Canada have a later phenology

# inspect random effects
re <- ranef(m_pheno_lat)$cond

levels <- list(
  list(label = "hex_id_res3",
    vars = c("hex_id_res3"),
    re = "hex_id_res3"),
  list(label = "hex_id_res4",
    vars = c("hex_id_res4","hex_id_res3"),
    re = "hex_id_res4:hex_id_res3"),
  list(label = "hex_id_res5",
    vars = c("hex_id_res5","hex_id_res4","hex_id_res3"),
    re = "hex_id_res5:hex_id_res4:hex_id_res3"),
  list(label = "hex_id_res6",
    vars = c("hex_id_res6","hex_id_res5","hex_id_res4","hex_id_res3"),
    re = "hex_id_res6:hex_id_res5:hex_id_res4:hex_id_res3"),
  list(label = "hex_id_res7",
    vars = c("hex_id_res7","hex_id_res6","hex_id_res5","hex_id_res4","hex_id_res3"),
    re = "hex_id_res7:hex_id_res6:hex_id_res5:hex_id_res4:hex_id_res3")
)

re_df <- map_dfr(levels, function(x){
  re_g <- re[[x$re]] %>%
    rownames_to_column("ID") %>%
    rename(intercept = `(Intercept)`)
  
  counts <- data_model %>%
    mutate(ID = do.call(paste, c(across(all_of(x$vars)), sep=":"))) %>%
    count(ID, name = "n_obs")
  
  left_join(re_g, counts, by = "ID") %>%
    mutate(level = x$label)
})

re_df <- re_df %>%
  mutate(level = factor(level, levels = c("hex_id_res3","hex_id_res4",
        "hex_id_res5","hex_id_res6","hex_id_res7"),
      labels = c("12,500 km2","750 km2","250 km2","90 km2", "25 km2")))

num_obs_vs_intercept <- ggplot(re_df, aes(n_obs, intercept)) +
  geom_jitter(
    width = 0.05,
    height = 0,
    alpha = 0.35,
    size = 1) +
  geom_smooth(method = "lm",formula = y ~ x,
    se = TRUE,
    color = "black",
    linewidth = 0.8) +
  geom_hline(yintercept = 0,
    linetype = "dashed",
    color = "red",
    linewidth = 0.6) +
  facet_wrap(~level, scales = "free") +
  scale_x_log10() +
  theme_bw() +
  labs(x = "Number of observations",
    y = "Estimated random intercept")

plot(num_obs_vs_intercept)

ggsave("diagnostics_and_checks/Plot - Number of observations vs. Random Intercept Estimate.png",
       plot   = num_obs_vs_intercept,
       width  = 30,
       height = 20,
       dpi    = 300,
       bg     = "white")

# no support for the hypothesis that sites with many observations were predicted
# as hotspots. Instead, hexbin random intercepts tended to decline when the hexagon
# contained large number of points
 
###
### VARIATION EXPLAINED BY THE MODEL COMPONENTS
###

library(performance)
library(insight)
library(glmmTMB)
library(tidyverse)
library(ggthemes)

# R2 calculations
r2 <- r2_nakagawa(m_pheno_lat)
r2_marginal <- r2$R2_marginal
r2_cond     <- r2$R2_conditional

# random effect variance
vc <- VarCorr(m_pheno_lat)$cond

var_hex7      <- vc$`hex_id_res7:hex_id_res6:hex_id_res5:hex_id_res4:hex_id_res3`[1,1]
var_hex6      <- vc$`hex_id_res6:hex_id_res5:hex_id_res4:hex_id_res3`[1,1]
var_hex5      <- vc$`hex_id_res5:hex_id_res4:hex_id_res3`[1,1]
var_hex4      <- vc$`hex_id_res4:hex_id_res3`[1,1]
var_hex3      <- vc$`hex_id_res3`[1,1]
var_year      <- vc$`YearCollected`[1,1]
var_ecoregion <- vc$`Ecoregion`[1,1]
var_survey    <- vc$`SurveyID`[1,1]

vars <- c(
  hex7      = var_hex7,
  hex6      = var_hex6,
  hex5      = var_hex5,
  hex4      = var_hex4,
  hex3      = var_hex3,
  year      = var_year,
  ecoregion = var_ecoregion,
  survey    = var_survey
)

# partitioning random variance
random_total <- r2_cond - r2_marginal
props <- vars / sum(vars)

parts_stack <- tibble(
  component = c(
    "Date*Latitude",
    "Year",
    "Ecoregion",
    "Survey",
    "Hex 12,500 km²",
    "Hex 750 km²",
    "Hex 250 km²",
    "Hex 90 km²",
    "Hex 25 km²",
    "Unexplained"),
  variance = c(
    r2_marginal,
    random_total * props["year"],
    random_total * props["ecoregion"],
    random_total * props["survey"],
    random_total * props["hex3"],
    random_total * props["hex4"],
    random_total * props["hex5"],
    random_total * props["hex6"],
    random_total * props["hex7"],
    1 - r2_cond))

parts_stack$component <- factor(parts_stack$component,
                                levels = rev(parts_stack$component))

# plot
cols <- c(
  "Hex 25 km²"     = "#67000D",
  "Hex 90 km²"     = "#A50F15",
  "Hex 250 km²"    = "#CB181D",
  "Hex 750 km²"    = "#EF3B2C",
  "Hex 12,500 km²" = "#FB6A4A",
  "Ecoregion"      = "darkblue",
  "Year"           = "gold",
  "Survey"         = "#8B4513",
  "Date*Latitude"  = "forestgreen",
  "Unexplained"    = "grey95")

parts_stack_plot <- ggplot(parts_stack,
                           aes(x = "Variance", y = variance, fill = component)) +
  geom_col(width = 0.8) +
  scale_fill_manual(values = cols, guide = guide_legend(reverse = FALSE)) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    expand = c(0, 0)) +
  labs(x = NULL, y = "Proportion of variance explained", fill = NULL) +
  theme_few(base_size = 13) +
  theme(
    axis.text.x        = element_blank(),
    axis.ticks.x       = element_blank(),
    panel.grid.major.x = element_blank(),
    aspect.ratio     = 3)

parts_stack_plot

###
### MAPPING OF HOTSPOTS
###

# random effects
re <- ranef(m_pheno_lat)$cond

# peak Julian day (at mean latitude) 
lat_mean  <- mean(data_model$latitude, na.rm = TRUE)

jday_grid <- data.frame(
  JulianDay     = seq(min(data_model$JulianDay), max(data_model$JulianDay), by = 1),
  latitude      = lat_mean,
  StemCount     = 1,
  SurveyID      = NA,                                          # added
  hex_id_res3   = NA, hex_id_res4 = NA, hex_id_res5 = NA,
  hex_id_res6   = NA, hex_id_res7 = NA,
  YearCollected = NA, Ecoregion   = NA
)

jday_grid$pred <- predict(m_pheno_lat, newdata = jday_grid,
                          type = "link", re.form = NA,
                          allow.new.levels = TRUE)

peak_jday <- jday_grid$JulianDay[which.max(jday_grid$pred)]

# fixed-effect prediction at peak Julian day
fe_at_peak_lat <- function(lat_vec) {
  grid <- data.frame(
    JulianDay     = peak_jday,
    latitude      = lat_vec,
    StemCount     = 1,
    SurveyID      = NA,                                        # added
    hex_id_res3   = NA, hex_id_res4 = NA, hex_id_res5 = NA,
    hex_id_res6   = NA, hex_id_res7 = NA,
    YearCollected = NA, Ecoregion   = NA
  )
  predict(m_pheno_lat, newdata = grid, type = "link",
          re.form = NA, allow.new.levels = TRUE)
}

# add centroid latitude to each hexgrid
add_latitude <- function(hex_sf) {
  hex_sf %>%
    mutate(latitude = st_coordinates(st_centroid(geometry))[, 2])
}

hex_res7 <- add_latitude(hex_res7)
hex_res6 <- add_latitude(hex_res6)
hex_res5 <- add_latitude(hex_res5)
hex_res4 <- add_latitude(hex_res4)
hex_res3 <- add_latitude(hex_res3)

# build per-resolution ranef data
make_re_sf <- function(re_component, res_level, hex_sf) {
  df <- data.frame(
    h3_index = sapply(strsplit(rownames(re_component), ":"), `[`, res_level),
    re_value = re_component[, 1]
  ) %>%
    group_by(h3_index) %>%
    summarise(re_value = mean(re_value))
  left_join(hex_sf, df, by = c("h3_address" = "h3_index"))
}

re_res7_sf <- make_re_sf(re[["hex_id_res7:hex_id_res6:hex_id_res5:hex_id_res4:hex_id_res3"]], 1, hex_res7)
re_res6_sf <- make_re_sf(re[["hex_id_res6:hex_id_res5:hex_id_res4:hex_id_res3"]],             1, hex_res6)
re_res5_sf <- make_re_sf(re[["hex_id_res5:hex_id_res4:hex_id_res3"]],                         1, hex_res5)
re_res4_sf <- make_re_sf(re[["hex_id_res4:hex_id_res3"]],                                     1, hex_res4)
re_res3_sf <- make_re_sf(re[["hex_id_res3"]],                                                 1, hex_res3)

# ecoregion ranef
re_eco_df <- data.frame(
  Ecoregion = rownames(re[["Ecoregion"]]),
  re_eco    = re[["Ecoregion"]][, 1])

hex_to_eco <- data_model %>%
  select(hex_id_res7, hex_id_res6, hex_id_res5, hex_id_res4, hex_id_res3, Ecoregion) %>%
  distinct()

# ecoregion RE and compute count_100 using per-hex latitude
add_eco_re <- function(sf_obj, hex_col) {
  sf_obj %>%
    left_join(
      hex_to_eco %>% select(all_of(hex_col), Ecoregion) %>% distinct(),
      by = c("h3_address" = hex_col)
    ) %>%
    left_join(re_eco_df, by = "Ecoregion") %>%
    mutate(
      re_value  = replace_na(re_value, 0),
      re_eco    = replace_na(re_eco, 0),
      fe_peak   = fe_at_peak_lat(latitude),
      count_100 = exp(fe_peak + re_value + re_eco + log(100))
    )
}

re_res7_sf <- add_eco_re(re_res7_sf, "hex_id_res7")
re_res6_sf <- add_eco_re(re_res6_sf, "hex_id_res6")
re_res5_sf <- add_eco_re(re_res5_sf, "hex_id_res5")
re_res4_sf <- add_eco_re(re_res4_sf, "hex_id_res4")
re_res3_sf <- add_eco_re(re_res3_sf, "hex_id_res3")

# midpoint reference
midpoint_count <- exp(fe_at_peak_lat(lat_mean) + log(100))

# plot over North America
library(rnaturalearth)
library(rnaturalearthdata)

north_america <- ne_countries(continent = "North America", returnclass = "sf") %>%
  st_transform(st_crs(hex_res7))

make_map3 <- function(sf_obj, title) {
  panel_range  <- range(log10(sf_obj$count_100), na.rm = TRUE)
  panel_breaks <- seq(panel_range[1], panel_range[2], length.out = 5)
  panel_labels <- round(10^panel_breaks, 1)

  ggplot() +
    geom_sf(data = north_america, fill = "black", colour = "grey60", linewidth = 0.2) +
    geom_sf(data = sf_obj, aes(fill = log10(count_100)), colour = NA) +
    scale_fill_gradient2(
      low      = "cyan", mid = "yellow", high = "red",
      midpoint = log10(midpoint_count),
      limits   = panel_range,
      breaks   = panel_breaks,
      labels   = panel_labels,
      name     = "Larvae per\n100 stems"
    ) +
    coord_sf(crs  = st_crs(hex_res7),
             xlim = st_bbox(sf_obj)[c("xmin", "xmax")],
             ylim = st_bbox(sf_obj)[c("ymin", "ymax")]) +
    labs(title = title) +
    theme_minimal() +
    theme(axis.text = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank())
}

make_map_points <- function(sf_obj, title) {
  panel_range  <- range(log10(sf_obj$count_100), na.rm = TRUE)
  panel_breaks <- seq(panel_range[1], panel_range[2], length.out = 5)
  panel_labels <- round(10^panel_breaks, 1)
  
  centroids <- sf_obj %>%
    mutate(geometry = st_centroid(geometry))
  
  ggplot() +
    geom_sf(data = north_america, fill = "black", colour = "grey60", linewidth = 0.2) +
    geom_sf(data = centroids, aes(colour = log10(count_100)), size = 0.3, shape = 16) +
    scale_colour_gradient2(
      low      = "cyan", mid = "yellow", high = "red",
      midpoint = log10(midpoint_count),
      limits   = panel_range,
      breaks   = panel_breaks,
      labels   = panel_labels,
      name     = "Larvae per\n100 stems"
    ) +
    coord_sf(crs  = st_crs(sf_obj),
             xlim = st_bbox(sf_obj)[c("xmin", "xmax")],
             ylim = st_bbox(sf_obj)[c("ymin", "ymax")]) +
    labs(title = title) +
    theme_minimal() +
    theme(axis.text = element_blank(), 
          axis.ticks = element_blank(), 
          axis.title = element_blank(),
          panel.grid = element_blank())+
    xlab(NULL) +
    ylab(NULL)
}

p7 <- make_map_points(re_res7_sf, "Hexagon area: ~25 km²")
p6 <- make_map_points(re_res6_sf, "Hexagon area: ~90 km²")
p5 <- make_map3(re_res5_sf, "Hexagon area: ~250 km²")
p4 <- make_map3(re_res4_sf, "Hexagon area: ~750 km²")
p3 <- make_map3(re_res3_sf, "Hexagon area: ~12,500 km²")

no_axes <- theme(
  axis.title = element_blank(),
  axis.text  = element_blank(),
  axis.ticks = element_blank()
)

p3 <- p3 + no_axes
p4 <- p4 + no_axes
p5 <- p5 + no_axes
p6 <- p6 + no_axes
p7 <- p7 + no_axes

pred_df <- data_model %>%
  mutate(pred_link  = predict(m_pheno_lat, type = "link",
                         re.form = NULL,
                         allow.new.levels = TRUE),
    pred_count = exp(pred_link) / StemCount * 100)

phist <- ggplot(pred_df, aes(x = pred_count)) +
  geom_histogram(bins = 50, fill = "orange", colour = "orange", linewidth = 0.2) +
  scale_x_log10(
    breaks = c(0.01, 0.1, 1, 10, 100, 1000),
    labels = c("0.01", "0.1", "1", "10", "100", "1000")
  ) +
  labs(x = "Predicted Monarch butterfly larvae per 100 stems", y = "Number of survey locations") +
  theme_minimal(base_size = 14) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"))


# final_plot <- (p7 + p6 + p5) / (p4 + p3 + phist) +
#   plot_annotation(
#     title = "Estimated larval counts per 100 stems\n(at peak phenology, hex-specific latitude)",
#     theme = theme(plot.title = element_text(face = "bold", size = 24))
#   )
# 
# final_plot

###
### YEAR RANDOM EFFECTS VS YEARLY TEMPERATURE
###

#since temperature is not currently in m_pheno_lat, this plot is asking whether unexplained year-to-year variation aligns with annual temperature.

# Extract year random effects
year_re_df <- data.frame(YearCollected = rownames(re[["YearCollected"]]),
  year_re = re[["YearCollected"]][, 1]) %>%
  mutate(YearCollected = as.factor(YearCollected))

# Calculate mean yearly temperature from the original model data
year_temp_df <- data %>%
  group_by(YearCollected) %>%
  summarise(mean_temp = mean(tmax_0y, na.rm = TRUE),   # change tmax to your temperature column name
    .groups = "drop")

# Join random effects to temperature
year_re_temp <- year_re_df %>%
  left_join(year_temp_df, by = "YearCollected")

# Plot relationship
p_year_temp <- ggplot(year_re_temp, aes(x = mean_temp, y = year_re)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = TRUE) +
  labs(x = "Mean yearly maximum temperature", y = "Year random effect") +
  theme_minimal(base_size = 14) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# alternative version
p_year_temp <- ggplot(year_re_temp, aes(x = mean_temp, y = year_re)) +
  geom_smooth(method = "lm", se = TRUE,
              colour = "darkorange", fill = "orange") +
  geom_point(size = 5, colour = "red") +
  labs(x = "Mean yearly maximum temperature", y = "Year random effect") +
  theme_minimal()+#)base_size = 54)+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    aspect.ratio     = 1.2   # height/width — increase to make it taller, decrease for wider
  )

# placeholder for monarch photo
placeholder <- ggplot() + theme_void() + theme(plot.background = element_rect(fill = "white", colour = NA))

# final_plot <- (parts_stack_plot + phist + placeholder) /
#   (p7 + p6 + p5) / 
#   (p4 + p3 + p_year_temp)  
# #final_plot
# 
# final_plot <- ((parts_stack_plot + phist + placeholder) /
#                  (p7 + p6 + p5) / 
#                  (p4 + p3 + p_year_temp)) +
#   plot_annotation(tag_levels = "A") &
#   theme(text         = element_text(size = 28),
#         axis.text    = element_text(size = 28),
#         legend.text  = element_text(size = 28),
#         title        = element_text(size = 28),
#         legend.title  = element_text(margin = margin(b = 24)),
#         legend.key.height = unit(24, "pt"),
#         legend.key.spacing  = unit(24, "pt")
#         )

library(patchwork)

design <- "
BCDE
AAAA
AAAA
AAAA
AAAA
FGGI
"
design <- "
AAAAB
AAAAC
AAAAD
AAAAE
FGGIH
"

final_plot <- (p3 + p4 + p5 + p6 + p7 +
                 parts_stack_plot + phist  + p_year_temp + placeholder) +
  plot_layout(design = design) +
  plot_annotation(tag_levels = "A") &
  theme(text         = element_text(size = 20),
        axis.text    = element_text(size = 20),
        legend.text  = element_text(size = 20),
        title        = element_text(size = 28),
        legend.title = element_text(margin = margin(b = 18)),
        legend.key.height  = unit(18, "pt"),
        legend.key.spacing = unit(18, "pt") )


ggsave("FIG2.png",
       plot   = final_plot,
       width  = 30,
       height = 20,
       dpi    = 300,
       bg     = "white")


