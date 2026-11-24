Code and data used to produce the results and figures found in Lewthwaite et. al. "Continental breeding hotspots of migratory monarch butterflies revealed by coordinated community science partnership"

This repository contains the scripts to complete all analyses in the manuscript, as well model outputs and required data inputs. Scripts are organized by the figure they correspond to, under the folder "scripts", and all data for both figures is in the folder "data". The model_outputs folder contains fitted statistical models, model comparison and diagnostic results, variance partitioning estimates, and spatial predictions used to generate Figure 2. These outputs allow users to inspect the model results and reproduce the figures without refitting the models.

# Figure 1: Descriptive Statistics for the Blitz Program

## Overview
This R script generates a three-panel figure summarizing the spatial and temporal characteristics of the Blitz monitoring program.

- **Panel A:** Geographic distribution of Blitz observations across North America.
- **Panel B:** Number of surveys and participants through time.
- **Panel C:** Temporal distribution of observations from Mission Monarch (MM) and the Monarch Larva Monitoring Project (MLMP), relative to Blitz onset (±30 days).

## Required Inputs
The script requires the following files:
1. "data/stats_long.csv" <- Annual numbers of surveys and participants.
2. "data/combined.csv" <- Observation data containing source, relative day, and proportion of observations.
3. "data/Mapping_points.csv" <- Geographic coordinates of Blitz observations.
4. "data/cat_wbg.png" <- Image incorporated into Panel C.

# Figure 2: Spatial and Temporal Variation in Monarch Larval Abundance

## Overview
This R script fits generalized linear mixed-effects models (GLMMs) to examine spatial and temporal variation in monarch larval abundance across North America. Larval counts are modeled using a negative binomial distribution, with Julian day, latitude, and their interaction as fixed effects, and survey, year, ecoregion, and five nested spatial resolutions as random effects. The number of milkweed stems surveyed is included as an offset.

The script compares alternative model distributions, evaluates model diagnostics, conducts sensitivity analyses, partitions variance among model components, and generates a multi-panel figure.

- **Panel A-E:**  Predicted monarch larval abundance per 100 milkweed stems across North America at five spatial resolutions (H3 hexagonal grids).
- **Panel F:**  Proportion of model-scale variance attributed to fixed effects, spatial random effects, year, ecoregion, survey, and unexplained variation.
- **Panel G:**  Distribution of predicted monarch larval abundance per 100 milkweed stems across survey locations.
- **Panel H:** Relationship between annual maximum temperature and year-specific random effects.

## Required Inputs

The script requires the following files:
1. "data/breeding_df.csv" <— Monarch breeding survey data, including larval counts, milkweed stem counts, observation dates, geographic coordinates, survey identifiers, ecoregions, and temperature variables.
2. "data/grid_50km/grid_50km.shp" — Spatial grid shapefile used to define the projected coordinate reference system. Associated shapefile files (e.g., .shx, .dbf, .prj) are also required.

H3 hexagonal grids are generated within the script and do not require separate input files.
