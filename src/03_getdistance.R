


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

# dist_df = as.data.frame(distance_mat) %>%
#         mutate(PARTICIPAN = rownames(distance_mat)) %>%
#         reshape2::melt(, id.vars = "PARTICIPAN") %>%
#         rename(other_pha = variable, 
#                distance = value) %>%
#         filter(!PARTICIPAN == other_pha) %>%
#         arrange(PARTICIPAN, distance) %>%
#         left_join(pha_nost %>% dplyr::select(PARTICIPAN, FORMAL_PAR) %>%
#                         rename(name_focal_pha = FORMAL_PAR),
#                   by = "PARTICIPAN") %>%
#         left_join(pha_nost %>% dplyr::select(PARTICIPAN, FORMAL_PAR) %>%
#                         rename(name_other_pha = FORMAL_PAR),
#                   by = c("other_pha" = "PARTICIPAN")) %>%
#         rename(geometry_focal = geometry.x,
#                geometry_other = geometry.y)
# 
# 
# fwrite(dist_df, "/Users/raj2/Dropbox/EligibilityPaper/PreferenceInference/Data/Raw/spatial_derivedobjs/pha_distance.csv")
# print("wrote file")
