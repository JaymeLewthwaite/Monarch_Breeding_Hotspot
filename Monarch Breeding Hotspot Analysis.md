Code and data used to produce the results and figures found in Lewthwaite et. al. "Continental breeding hotspots of migratory monarch butterflies revealed by coordinated community science partnership"

This repository contains the scripts to complete all analyses in the manuscript, as well model outputs (???) and required data inputs are included. Scripts are organized by the figure they correspond to, under the folder "scripts", and all data for both figures is in the folder "data"

#Figure 1: Descriptive Statistics for the Blitz Program

##Overview
This R script generates a three-panel figure summarizing the spatial and temporal characteristics of the Blitz monitoring program.

Panel A: Geographic distribution of Blitz observations across North America.
Panel B: Number of surveys and participants through time.
Panel C: Temporal distribution of observations from Mission Monarch (MM) and the Monarch Larva Monitoring Project (MLMP), relative to Blitz onset (±30 days).

Required Inputs
The script requires the following files:

1. "data/stats_long.csv" <- Annual numbers of surveys and participants.

2. "data/combined.csv" <- Observation data containing source, relative day, and proportion of observations.

3. "data/Mapping_points.csv" <- Geographic coordinates of Blitz observations.

4. "cat_wbg.png" <- Image incorporated into Panel C.
