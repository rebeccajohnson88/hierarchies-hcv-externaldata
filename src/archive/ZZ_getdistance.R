


## load packages
library(tigris)
library(ggplot2)
library(rgdal)
library(sf)
library(tigris)
library(tidyverse)
library(tictoc)

## read in pha shapefile
options(tigris_class = "sf")
options(tigris_use_cache = TRUE)


pha_df = readRDS("/Users/raj2/Dropbox/EligibilityPaper/PreferenceInference/Data/Raw/spatial_derivedobjs/phas_foroverlap.RDS")

## exclude ones that are statewide
pha_nost = pha_df %>%
      mutate(statewide =  ifelse(grepl("^[A-Z][A-Z]9",
               PARTICIPAN), 1, 0)) %>%
      filter(statewide == 0)

## get centroids (makes distance calc faster)
pha_nost_cent = st_centroid(pha_nost)

## start clock
tic("distance calc")

## test with small sample
#pha_nost_cent_sample = pha_nost_cent %>%
 #               filter(STUSAB == "AL")
print("starting distance calc")
distance_mat = st_distance(pha_nost_cent, pha_nost_cent)
rownames(distance_mat) = pha_nost_cent$PARTICIPAN
colnames(distance_mat) = pha_nost_cent$PARTICIPAN
print("ended distance calc")
saveRDS(distance_mat, 
    "/Users/raj2/Dropbox/EligibilityPaper/PreferenceInference/Data/Raw/spatial_derivedobjs/centroid_distance_mat.RDS")

