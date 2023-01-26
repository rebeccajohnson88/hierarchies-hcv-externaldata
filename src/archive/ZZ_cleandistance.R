


## load packages
library(tigris)
library(ggplot2)
library(rgdal)
library(sf)
library(tigris)
library(tidyverse)
library(tictoc)
library(data.table)

## read in pha shapefile
options(tigris_class = "sf")
options(tigris_use_cache = TRUE)

distance_mat = readRDS("/Users/raj2/Dropbox/EligibilityPaper/PreferenceInference/Data/Raw/spatial_derivedobjs/centroid_distance_mat.RDS")
pha_df = readRDS("/Users/raj2/Dropbox/EligibilityPaper/PreferenceInference/Data/Raw/spatial_derivedobjs/phas_foroverlap.RDS")

## exclude ones that are statewide
pha_nost = pha_df %>%
  mutate(statewide =  ifelse(grepl("^[A-Z][A-Z]9",
                                   PARTICIPAN), 1, 0)) %>%
  filter(statewide == 0)


dist_df = as.data.frame(distance_mat) %>%
         mutate(PARTICIPAN = rownames(distance_mat)) %>%
        reshape2::melt(, id.vars = "PARTICIPAN") %>%
         rename(other_pha = variable, 
               distance = value) %>%
         filter(!PARTICIPAN == other_pha) %>%
        arrange(PARTICIPAN, distance) 

fwrite(dist_df, "/Users/raj2/Dropbox/EligibilityPaper/PreferenceInference/Data/Raw/spatial_derivedobjs/pha_distance.csv")

