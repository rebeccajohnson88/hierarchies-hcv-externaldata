# 00a_pulltract_shapefiles.R
# Script that uses the 
# sf package to pull tract-level shape files
# for each of the states in the raw PHA data

library(here)
source(here("src/utils.R"))

WRITE_NEW_INTERMEDIATE_FILES <- FALSE 

###########################################################################
## Get stusab codes from raw service area data
###########################################################################

pha_shpformat <- readOGR(here("data/raw/Estimated_Housing_Authority_Service_Areas/"))
state_codes = unique(pha_shpformat$STUSAB)

###########################################################################
## Get tract shapefiles for each of the states
###########################################################################

tracts_allstates = lapply(state_codes,
                function(x) tracts(x, cb = TRUE))
tracts_merged = rbind_tigris(tracts_allstates)

###########################################################################
## Save in RDS format
###########################################################################

if(WRITE_NEW_INTERMEDIATE_FILES){
  saveRDS(tracts_merged,
          here("data/intermediate/tracts_sf_format.RDS"))
}


