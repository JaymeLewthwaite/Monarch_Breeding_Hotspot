# Figure 1. Descriptive statistics for the Blitz program

# ======================================================================
# Packages
# ======================================================================

library(tidyverse)
library(sf)
library(rnaturalearth)
library(png)

# import data
stats_long <- read_csv("data/stats_long.csv")
combined    <- read_csv("data/combined.csv")
points <- read_csv("data/Mapping_points.csv")

# ======================================================================
# Prepare data for Panel A
# ======================================================================

points_sf <- st_as_sf(
  points,
  coords = c("DecimalLongitude", "DecimalLatitude"),
  crs = 4326
)

north_america <- ne_countries(
  scale = "medium",
  returnclass = "sf"
) #%>%
  #filter(continent == "North America")

bbox <- st_bbox(points_sf)

x_margin <- 0.01 * (bbox["xmax"] - bbox["xmin"])
y_margin <- 0.01 * (bbox["ymax"] - bbox["ymin"])


# ======================================================================
# Prepare data for Panel B
# ======================================================================

n_survey <- stats_long %>%
  filter(variable == "Number of surveys")

n_participant <- stats_long %>%
  filter(variable == "Number of participants")

# ======================================================================
# Prepare data for Panel C
# ======================================================================

combined_30 <- combined %>%
  filter(relative_day >= -30,
         relative_day <= 30)

mm_data <- combined_30 %>%
  filter(source == "MM")

mlmp_data <- combined_30 %>%
  filter(source == "MLMP")

# Function for calculating min/max intervals
get_range <- function(df) {
  tmp <- tapply(df$pct,
                factor(df$relative_day, levels = -30:30),
                range,
                na.rm = TRUE)
  mat <- matrix(unlist(tmp), ncol = 2, byrow = TRUE)
  t(apply(mat, 1, sort))
}

mm_range   <- get_range(mm_data)
mlmp_range <- get_range(mlmp_data)

# Colors
col_mm   <- rgb(219, 193, 42, maxColorValue = 255, alpha = 150)
col_mlmp <- rgb(27, 171, 63, maxColorValue = 255, alpha = 150)


# ======================================================================
# Figure
# ======================================================================
#pdf("figure1.pdf", width = 9, height = 3)

png("figure1.png", width = 9, height = 3.5, units = "in", res = 300)

layout(
  matrix(c(1, 2, 3), nrow = 1),
  widths = c(1.15, 1, 1)   # Panel A plus large
)

par(
  pty = "m",
  oma = c(0, 0, 0, 0),
  mgp = c(2, 0.6, 0),
  tcl = -0.25
)

# ----------------------------------------------------------------------
# Panel A. Spatial distribution of Blitz observations
# ----------------------------------------------------------------------
par(mar = c(4, 1, 6, 1))

plot(
  st_geometry(north_america),
  col = "black",
  border = "black",
  xlim = c(
    bbox["xmin"],
    bbox["xmax"]
  ),
  ylim = c(
    bbox["ymin"],
    bbox["ymax"]
  ),
  axes = FALSE
)

box(col = "black")

# Points
plot(st_geometry(points_sf),
     add = TRUE,
     col = "#fc8d62",
     pch = 16,
     cex = 0.3)

mtext("A", 3, adj = 0, cex = 2)


# ----------------------------------------------------------------------
# Panel B. Number of surveys and participants through time
# ----------------------------------------------------------------------
par(mar = c(4, 3, 6, 0))

plot(
  n_survey$YearCollected,
  n_survey$value,
  type = "l",
  col = "#fc8d62",
  lty = 1,
  lwd = 3,
  pch = 19,
  xaxt = "n",
  yaxt = "n",
  xlab = "",
  ylab = ""
)

axis(1, at = nsurvey$YearCollected, las = 2)
axis(2, las = 1)

lines(
  n_participant$YearCollected,
  n_participant$value,
  type = "l",
  col = "#fc8d62",
  lty = 3,
  lwd = 3,
  pch = 19)

legend("topleft",
       inset = c(0.05, 0),
       legend = c("Number of surveys", "Number of participants"),
       col = "#fc8d62",
       pch = NA,
       lty = c(1, 3),
       lwd = 3,
       bty = "n"
)

mtext("B", 3, adj = 0, cex = 2)

#----------------
# PANEL C
#----------------
par(mar = c(4, 4, 6, 1), mgp = c(2.5, 0.6, 0), xpd = NA)

plot(
  NULL,
  xlim = c(-30, 30),
  ylim = c(
    0,
    max(mm_range, mlmp_range, na.rm = TRUE) * 1.05
  ),
  xlab = "Days relative to Blitz onset",
  ylab = "Proportion of observations",
  las = 1
)

xvals <- -30:30

polygon(
  c(xvals, rev(xvals)),
  c(mm_range[, 2], rev(mm_range[, 1])),
  border = NA,
  col = col_mm
)

polygon(
  c(xvals, rev(xvals)),
  c(mlmp_range[, 2], rev(mlmp_range[, 1])),
  border = NA,
  col = col_mlmp
)

abline(
  v = c(0, 10),
  col = adjustcolor("darkgrey", alpha.f = 0.4),
  lwd = 3,
  xpd = FALSE
)

legend(
  "topleft",
  legend = c("MM", "MLMP"),
  fill = c(col_mm, col_mlmp),
  bty = "n"
)

mtext("C", 3, adj = 0, cex = 2)

img <- readPNG("cat_wbg.png")

usr <- par("usr")

rasterImage(
  img,
  xleft   = -20,
  ybottom = usr[4] * 1.003,
  xright  = 32.55,
  ytop    = usr[4] * 1.3
)

dev.off()